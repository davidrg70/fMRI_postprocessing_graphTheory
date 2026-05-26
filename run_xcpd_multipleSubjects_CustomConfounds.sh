#!/usr/bin/env bash
set -euo pipefail

module purge
module load apptainer

# ---------------- Paths ----------------
XCP_DIR="/work/users/d/g/dga/tools/XCP-D"
SIMG="${XCP_DIR}/XCP-D-0.10.7.simg"

FMRIPREP_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep"
OUTPUT_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd"
ATLAS_DIR="/proj/cohenlab/projects/ADHDBrainMAP/BrainMap_Proc_clpipe191_Nov2024/Seitzman"
CUSTOM_DATASET="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep/custom_confounds"
NUISANCE_YAML="/users/d/g/dga/BrainMAP/xcpd_work/custom_config.yaml"

LOG_DIR="${OUTPUT_DIR}/slurm_logs"
mkdir -p "${LOG_DIR}"

# --- Participants to process ---
PARTICIPANT_LABELS=('IDS_IN_HERE')

# --- SLURM resources ---
CPUS=4
MEM="128G"
TIME="08:00:00"
PARTITION="small"   # on Sycamore partition

# --- Task family you care about ---
TASK_PREFIX="gngregular"   # will match gngregular01, gngregular02, etc.

for PART in "${PARTICIPANT_LABELS[@]}"; do
  FUNC_DIR="${FMRIPREP_DIR}/sub-${PART}/ses-01/func"

  # Find candidate BOLD files in standard space with desc-preproc
  mapfile -t BOLDS < <(ls "${FUNC_DIR}"/sub-"${PART}"*_space-MNI152NLin2009cAsym*_desc-preproc_bold.nii.gz 2>/dev/null || true)

  if [[ ${#BOLDS[@]} -eq 0 ]]; then
    echo "!! No MNI preproc BOLD found for ${PART} in ${FUNC_DIR} (skipping)"
    continue
  fi

  # Extract unique task labels from filenames (robust to run formatting)
  # Pull the substring after 'task-' until the next '_' character
  TASKS=$(printf "%s\n" "${BOLDS[@]}" | sed -n 's/.*_task-\([^_]*\)_.*/\1/p' | sort -u)

  # Keep only tasks that start with gngregular
  TASKS=$(printf "%s\n" "${TASKS}" | awk -v pfx="${TASK_PREFIX}" '$0 ~ ("^" pfx) {print}')

  if [[ -z "${TASKS}" ]]; then
    echo "!! No tasks matching ${TASK_PREFIX} for ${PART} (skipping)"
    echo "   Found tasks were:"
    printf "   - %s\n" $(printf "%s\n" "${BOLDS[@]}" | sed -n 's/.*_task-\([^_]*\)_.*/\1/p' | sort -u)
    continue
  fi

  while read -r TASK_ID; do
    [[ -z "${TASK_ID}" ]] && continue

    JOB_NAME="xcpd_${PART}_${TASK_ID}"
    JOB_SCRIPT="${LOG_DIR}/job_${JOB_NAME}.sh"

    cat > "${JOB_SCRIPT}" <<EOT
#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH -J ${JOB_NAME}
#SBATCH --output=${LOG_DIR}/${JOB_NAME}_%j.out
#SBATCH --error=${LOG_DIR}/${JOB_NAME}_%j.err
#SBATCH --mail-user=dga@ad.unc.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
module load apptainer

PART="${PART}"
TASK_ID="${TASK_ID}"

WORK_DIR="${OUTPUT_DIR}/work_\${PART}_\${TASK_ID}"
mkdir -p "\${WORK_DIR}" "${OUTPUT_DIR}"

NCPUS="\${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="\${NCPUS}"

BINDPATHS="/users,/proj,/work"

apptainer run --cleanenv -B "\${BINDPATHS}" \\
  "${SIMG}" \\
  "${FMRIPREP_DIR}" \\
  "${OUTPUT_DIR}" \\
  participant \\
  --fs-license-file "/work/users/d/g/dga/tools/freesurfer/7.3.2/license.txt" \\
  --mode nichart \\
  --participant-label "\${PART}" \\
  -d atlas="${ATLAS_DIR}" custom="${CUSTOM_DATASET}" \\
  --atlases Seitzman \\
  -p "${NUISANCE_YAML}" \\
  -t "\${TASK_ID}" \\
  --file-format nifti \\
  --smoothing 0 \\
  --motion-filter-type none \\
  -f 0.3 \\
  --min-time 60 \\
  --output-type interpolated \\
  --lower-bpf 0.01 \\
  --upper-bpf 0.08 \\
  --min-coverage 0.5 \\
  --nprocs "\${NCPUS}" \\
  --clean-workdir \\
  -w "\${WORK_DIR}"
EOT

    sed -i 's/\r$//' "${JOB_SCRIPT}"
    chmod +x "${JOB_SCRIPT}"

    echo "Submitting ${JOB_NAME}"
    sbatch "${JOB_SCRIPT}"

  done <<< "${TASKS}"
done

echo "✅ Done submitting."

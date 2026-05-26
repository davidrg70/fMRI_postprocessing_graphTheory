#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH -J xcpd_2002W01
#SBATCH --output=xcpd_idToChangeHere_%j.out
#SBATCH --error=xcpd_idToChangeHere_%j.err
#SBATCH --mail-user=dga@ad.unc.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load apptainer

# ---------------- Paths ----------------
XCP_DIR="/work/users/d/g/dga/tools/XCP-D"
SIMG="${XCP_DIR}/XCP-D-0.10.7.simg"

FMRIPREP_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep"
OUTPUT_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/"
ATLAS_DIR="/proj/cohenlab/projects/ADHDBrainMAP/BrainMap_Proc_clpipe191_Nov2024/Seitzman"
CUSTOM_DATASET="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep/custom_confounds"
NUISANCE_YAML="/users/d/g/dga/BrainMAP/xcpd_work/custom_config.yaml"

PART="2002W01"
TASK_ID="gngregular"

# Subject-specific work dir to avoid collisions
WORK_DIR="${OUTPUT_DIR}/work_${PART}"

mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

# Use allocated cores
NCPUS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${NCPUS}"

# Bind needed host paths into the container
BINDPATHS="/users,/proj,/work"

apptainer run --cleanenv -B "${BINDPATHS}" \
  "${SIMG}" \
  "${FMRIPREP_DIR}" \
  "${OUTPUT_DIR}" \
  participant \
  --fs-license-file "/work/users/d/g/dga/tools/freesurfer/7.3.2/license.txt" \
  --mode nichart \
  --participant-label "${PART}" \
  -d atlas="${ATLAS_DIR}" custom="${CUSTOM_DATASET}" \
  --atlases Seitzman \
  -p "${NUISANCE_YAML}" \
  -t "${TASK_ID}" \
  --file-format nifti \
  --smoothing 0 \
  --motion-filter-type none \
  -f 0.3 \
  --min-time 60 \
  --output-type interpolated \
  --lower-bpf 0.01 \
  --upper-bpf 0.08 \
  --min-coverage 0.5 \
  --nprocs "${NCPUS}" \
  --clean-workdir \
  -w "${WORK_DIR}"

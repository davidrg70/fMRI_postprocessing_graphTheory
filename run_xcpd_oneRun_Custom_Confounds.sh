#!/usr/bin/env bash
set -euo pipefail

module purge
module load apptainer

XCP_DIR="/work/users/d/g/dga/tools/XCP-D"
SIMG="${XCP_DIR}/XCP-D-0.10.7.simg"
FMRIPREP_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep/"
OUTPUT_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/"
CUSTOM_DATASET="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_fmriprep/custom_confounds"
NUISANCE_YAML="/users/d/g/dga/BrainMAP/xcpd_work/custom_config.yaml"
PART="2002W01"
TASK_ID="gngregular"
SESSION_ID="01"
WORK_DIR="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/work_2002W01"
BIDS_FILTER="/users/d/g/dga/BrainMAP/xcpd_work/data_gngregular_xcpd/run01_filter.json"

mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}"

NCPUS=2
export OMP_NUM_THREADS="${NCPUS}"
BINDPATHS="/users,/proj,/work"

apptainer run --cleanenv -B "${BINDPATHS}" \
  "${SIMG}" \
  "${FMRIPREP_DIR}" \
  "${OUTPUT_DIR}" \
  participant \
  --fs-license-file "/work/users/d/g/dga/tools/freesurfer/7.3.2/license.txt" \
  --mode nichart \
  --participant-label "${PART}" \
  --session-id "${SESSION_ID}" \
  --bids-filter-file "${BIDS_FILTER}" \
  -d custom="${CUSTOM_DATASET}" \
  -p "${NUISANCE_YAML}" \
  -t "${TASK_ID}" \
  --file-format nifti \
  --smoothing 0 \
  --motion-filter-type none \
  -f 0.3 \ # this will run censoring: scrubbing
  --min-time 60 \
  --output-type interpolated \ # this will run interpolation
  --lower-bpf 0.01 \
  --upper-bpf 0.08 \
  --min-coverage 0.5 \
  --nprocs "${NCPUS}" \
  --clean-workdir \
  -w "${WORK_DIR}"

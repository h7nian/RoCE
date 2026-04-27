#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=log/rhc_download_%j.out
#SBATCH --error=log/rhc_download_%j.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --job-name=FACEC_rhc_dl
#SBATCH -p msismall,amdsmall,agsmall
#SBATCH --nice=5

# ============================================================================
# FACE-C RHC Public Dataset Download
# ============================================================================
# Fetches the Right-Heart Catheterization teaching dataset (Connors et al.
# 1996) from Frank Harrell's hbiostat.org mirror into the package's
# inst/extdata/ directory.
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

SRC_URL="https://hbiostat.org/data/repo/rhc.csv"
DEST_PATH="inst/extdata/rhc.csv"

mkdir -p "$(dirname "${DEST_PATH}")" log

echo "=============================================="
echo "FACE-C RHC Download Job"
echo "=============================================="
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Source:      ${SRC_URL}"
echo "Destination: ${DEST_PATH}"
echo "Start time:  $(date)"
echo "=============================================="

# wget with:
#   -q quiet (all status goes to SLURM logs via --show-progress)
#   --show-progress makes the single-file progress visible
#   --tries=3 --timeout=60 for minor transient network failures
#   -O forces the output path even if intermediate redirects
wget --show-progress --tries=3 --timeout=60 -O "${DEST_PATH}" "${SRC_URL}"
WGET_EXIT=$?

if [[ ${WGET_EXIT} -ne 0 ]]; then
    echo "wget failed (exit ${WGET_EXIT})"
    exit ${WGET_EXIT}
fi

if [[ ! -s "${DEST_PATH}" ]]; then
    echo "Download completed but ${DEST_PATH} is empty"
    exit 2
fi

BYTES=$(stat -c '%s' "${DEST_PATH}")
LINES=$(wc -l < "${DEST_PATH}")

echo ""
echo "=============================================="
echo "Download complete."
echo "  File:  ${DEST_PATH}"
echo "  Size:  ${BYTES} bytes"
echo "  Lines: ${LINES} (including header)"
echo "  Time:  $(date)"
echo "=============================================="

#!/bin/bash
# ============================================================================
# Download the Right-Heart Catheterization (RHC) Public Dataset
# ============================================================================
# One-shot SLURM submission to fetch the RHC teaching dataset published by
# Frank Harrell at hbiostat.org and vendor it into the package as
# inst/extdata/rhc.csv. The dataset is ~1 MB (5,735 rows) and originates
# from Connors et al. (1996) NEJM 334:1441-1448.
#
# Usage:
#   ./data/download_rhc.sh                # Submit the SLURM download job
#   ./data/download_rhc.sh --dry-run      # Show sbatch command without submitting
#
# Once this script has been run successfully, inst/extdata/rhc.csv is
# available to any installation of the FACEHD package (via load_rhc_raw()
# and related helpers in R/real_data_rhc.R).
# ============================================================================

set -euo pipefail

DRY_RUN=false
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        -h|--help)
            sed -n '2,18p' "$0"; exit 1 ;;
        *) echo "Error: unknown option '$arg'"; exit 1 ;;
    esac
done

mkdir -p log

echo "======================================================"
echo " FACE-HD Real-Data Download: RHC"
echo "======================================================"
echo " Source:      https://hbiostat.org/data/repo/rhc.csv"
echo " Destination: inst/extdata/rhc.csv"
echo " Dry run:     ${DRY_RUN}"
echo "======================================================"

if [[ "$DRY_RUN" == "true" ]]; then
    echo ""
    echo "[DRY RUN] Would submit: sbatch --job-name=FACEHD_rhc_dl data/download_rhc.cmd"
    echo ""
    exit 0
fi

submit_output=$(sbatch --job-name=FACEHD_rhc_dl data/download_rhc.cmd)
echo ""
echo "Submitted: ${submit_output}"
echo ""
echo " Useful commands:"
echo "   squeue -u \$USER"
echo "   tail -f log/rhc_download_*.out"
echo "======================================================"

#!/bin/bash
# ============================================================================
# FACE-HD Real-Data (RHC) Job Submission Script
# ============================================================================
#
# Submits a single SLURM job that runs the Right-Heart Catheterization (RHC)
# real-data experiment (realdata.R / realdata.cmd). Parallel in spirit to
# main.sh for the simulation study, but without a parameter grid — the RHC
# experiment is a single-shot run parameterised only by K, outcome, n_folds
# and A_val.
#
# Usage:
#   ./realdata.sh                              # K=5, outcome=death30, n_folds=10, A_val=1
#   ./realdata.sh --outcome death180           # switch outcome to 180-day mortality
#   ./realdata.sh --outcome los                # continuous outcome (hospital LOS)
#   ./realdata.sh --K 4                        # 4 sites instead of the 5 default
#   ./realdata.sh --n-folds 5                  # lower-K cross-fitting
#   ./realdata.sh --A-val 0                    # estimate mu^0 arm instead of mu^1
#   ./realdata.sh --target-site CHF            # override default (largest cat1) target
#   ./realdata.sh --dry-run                    # print the sbatch command without submitting
#
# Directory layout:
#   log/                      SLURM stdout/stderr (realdata_<jobid>.out/.err)
#   results/real_data/        CSV tables + RDS + PDF figures
#
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# 1. Defaults
# ----------------------------------------------------------------------------
K=5
OUTCOME="death30"
N_FOLDS=10
A_VAL=1
SITE_VAR="cat1"
M_TAU_INFERENCE="3"
SITE_PRESET=""
PHI="identity"
TARGET_SITE=""
DRY_RUN=false

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --K N                   Number of sites (target + sources), default 5"
    echo "  --outcome TYPE          death30 (default) | death180 | los"
    echo "  --n-folds N             Outer cross-fitting folds, default 10"
    echo "  --A-val 0|1             Target treatment level, default 1 (RHC exposure)"
    echo "  --site-var NAME         Raw RHC column used for site partitioning"
    echo "                           (default 'cat1'; alternatives: ninsclas, income, race, ca)"
    echo "  --m-tau-inference VAL   Inference-time truncation cap on linear predictor"
    echo "                           (default 3; accepts 'Inf' for no truncation)"
    echo "  --site-preset NAME      Optional recode of site_var levels before partitioning."
    echo "                           Available presets:"
    echo "                             ninsclas4  — drop 'No insurance'; merge 'Medicaid'"
    echo "                                          + 'Medicare & Medicaid' → 'Medicaid_any'"
    echo "                                          (4 levels total; use with --site-var ninsclas)"
    echo "  --phi NAME              Covariate basis map applied inside build_rhc_data_split()."
    echo "                           Available values:"
    echo "                             identity  — default; no dimension inflation"
    echo "                             bspline   — cubic B-spline (knots=5) on continuous vars"
    echo "                                          (p rises from ~75 to ~200; ell_1 penalty active)"
    echo "  --target-site NAME      Override target site-level value (default: largest category)"
    echo "  --dry-run               Print sbatch command without submitting"
    echo "  -h, --help              Show this message"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --K)                K="${2:?Missing value for --K}"; shift 2 ;;
        --outcome)          OUTCOME="${2:?Missing value for --outcome}"; shift 2 ;;
        --n-folds)          N_FOLDS="${2:?Missing value for --n-folds}"; shift 2 ;;
        --A-val)            A_VAL="${2:?Missing value for --A-val}"; shift 2 ;;
        --site-var)         SITE_VAR="${2:?Missing value for --site-var}"; shift 2 ;;
        --m-tau-inference)  M_TAU_INFERENCE="${2:?Missing value for --m-tau-inference}"; shift 2 ;;
        --site-preset)      SITE_PRESET="${2:?Missing value for --site-preset}"; shift 2 ;;
        --phi)              PHI="${2:?Missing value for --phi}"; shift 2 ;;
        --target-site)      TARGET_SITE="${2:?Missing value for --target-site}"; shift 2 ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)          usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

# ----------------------------------------------------------------------------
# 2. Validate
# ----------------------------------------------------------------------------
if ! [[ "$K" =~ ^[0-9]+$ ]] || [[ "$K" -lt 2 ]]; then
    echo "Error: --K must be an integer >= 2 (got '$K')"; exit 1
fi
if ! [[ "$N_FOLDS" =~ ^[0-9]+$ ]] || [[ "$N_FOLDS" -lt 3 ]]; then
    echo "Error: --n-folds must be an integer >= 3 (got '$N_FOLDS')"; exit 1
fi
if [[ "$A_VAL" != "0" && "$A_VAL" != "1" ]]; then
    echo "Error: --A-val must be 0 or 1 (got '$A_VAL')"; exit 1
fi
case "$OUTCOME" in
    death30|death180|los) ;;
    *) echo "Error: --outcome must be one of death30, death180, los (got '$OUTCOME')"; exit 1 ;;
esac
case "$PHI" in
    identity|bspline) ;;
    *) echo "Error: --phi must be 'identity' or 'bspline' (got '$PHI')"; exit 1 ;;
esac

# ----------------------------------------------------------------------------
# 3. Summary
# ----------------------------------------------------------------------------
mkdir -p log results/real_data

# Build a setting_id suffix "Mt<n>" / "Mtinf" that keeps M_tau variants
# from overwriting each other's output files.
case "${M_TAU_INFERENCE}" in
    Inf|inf|Infinity|infinity) M_TAU_TAG="Mtinf" ;;
    *)                          M_TAU_TAG="Mt${M_TAU_INFERENCE/./p}" ;;  # e.g. 2.5 -> Mt2p5
esac

SETTING_ID="rhc_K${K}_${OUTCOME}_kf${N_FOLDS}_A${A_VAL}_${SITE_VAR}_${M_TAU_TAG}"
if [[ -n "$SITE_PRESET" ]]; then
    SETTING_ID="${SETTING_ID}_${SITE_PRESET}"
fi
if [[ "$PHI" != "identity" ]]; then
    SETTING_ID="${SETTING_ID}_${PHI}"
fi

echo "======================================================"
echo " FACE-HD Real-Data (RHC) Job Submission"
echo "======================================================"
echo " Setting ID:    ${SETTING_ID}"
echo " K:             ${K}"
echo " Outcome:       ${OUTCOME}"
echo " n_folds:       ${N_FOLDS}"
echo " A_val:         ${A_VAL}"
echo " Site var:      ${SITE_VAR}"
echo " M_tau_inf:     ${M_TAU_INFERENCE}"
echo " Site preset:   ${SITE_PRESET:-<none>}"
echo " Phi basis:     ${PHI}"
echo " Target site:   ${TARGET_SITE:-<auto: largest (possibly remapped) level of ${SITE_VAR}>}"
echo " Dry run:       ${DRY_RUN}"
echo "======================================================"

# ----------------------------------------------------------------------------
# 4. Submit (or print)
# ----------------------------------------------------------------------------
EXPORT_VARS="K_arg=${K},OUTCOME=${OUTCOME},N_FOLDS=${N_FOLDS},A_VAL=${A_VAL},SITE_VAR=${SITE_VAR},M_TAU_INFERENCE=${M_TAU_INFERENCE},SITE_PRESET=${SITE_PRESET},PHI=${PHI},TARGET_SITE=${TARGET_SITE}"

if [[ "$DRY_RUN" == "true" ]]; then
    echo ""
    echo "[DRY RUN] Would submit:"
    echo "  sbatch --export=${EXPORT_VARS} --job-name=FACE_realdata realdata.cmd"
    echo ""
else
    submit_output=$(sbatch \
        --export="${EXPORT_VARS}" \
        --job-name="FACE_realdata_${OUTCOME}_${M_TAU_TAG}" \
        realdata.cmd)
    echo ""
    echo "Submitted: ${submit_output}"
fi

echo ""
echo "======================================================"
if [[ "$DRY_RUN" == "true" ]]; then
    echo " DRY RUN COMPLETE — nothing submitted"
else
    echo " DONE — real-data job submitted"
fi
echo "======================================================"
echo ""
echo " Useful commands:"
echo "   squeue -u \$USER                          # monitor queue"
echo "   tail -f log/realdata_*.out                # watch live output"
echo "   ls results/real_data/                     # list outputs"
echo "======================================================"

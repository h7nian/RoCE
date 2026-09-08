#!/bin/bash
#SBATCH --job-name=refresh_figs
#SBATCH --time=00:25:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=2
#SBATCH --output=diagnosis/face_probe/validation/face_refresh_figs_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
echo "== re-aggregate negative-transfer chunks =="
Rscript diagnosis/face_probe/face_negT_aggregate.R 2>&1 | tail -3
echo "== plot p50 + p100 (C1/C2/C3) =="
for P in 50 100; do for cfg in C1 C2 C3; do
  export ROCE_FIG_CONFIG=$cfg ROCE_FIG_P=$P
  Rscript diagnosis/face_probe/face_fig1_plot.R 2>&1 | grep wrote
done; done
echo "== sync to docs/figures =="
for cfg in C1 C2 C3; do for P in 50 100; do
  cp diagnosis/face_probe/validation/fig1_negtransfer_${cfg}_p${P}.pdf docs/figures/sim_negtransfer_${cfg}_p${P}.pdf 2>/dev/null
done; done
echo "== continuous 2-round re-plot (p10+p50) =="
Rscript diagnosis/face_probe/face_2round_cont_plot.R 2>&1 | tail -2
echo "DONE"

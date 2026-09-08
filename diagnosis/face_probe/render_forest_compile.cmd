#!/bin/bash
#SBATCH --job-name=forest_compile
#SBATCH --time=00:25:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=2
#SBATCH --output=diagnosis/face_probe/validation/forest_compile_%j.out
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS=$HOME/Rlibs_em:$HOME/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
echo "===== rendering two-panel RHC forest (mu^1 + mu^0) ====="
Rscript diagnosis/face_probe/forest_replot_both.R
echo "render exit: $?"
ls -la docs/figures/rhc_forest_both.pdf 2>/dev/null || echo "[FAIL] forest pdf not written"

module load texlive/2025 2>/dev/null
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"/docs
for f in main supplemental; do
  echo "===== compiling $f.tex ====="
  pdflatex -interaction=nonstopmode -halt-on-error $f.tex >/tmp/${f}_c1.log 2>&1
  bibtex $f >/tmp/${f}_bib.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error $f.tex >/tmp/${f}_c2.log 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error $f.tex >/tmp/${f}_c3.log 2>&1
  if [ -f $f.pdf ]; then
    echo "[OK] $f.pdf  ($(pdfinfo $f.pdf 2>/dev/null | grep -i Pages | tr -s ' '))"
  else
    echo "[FAIL] $f.tex did not produce a PDF"; grep -iE '^!|error' /tmp/${f}_c1.log | head -10
  fi
  echo "--- errors (last pass) ---"
  grep -iE '^!|! LaTeX Error|Undefined control|Runaway|Missing' /tmp/${f}_c3.log | head -8
  echo "--- undefined refs/citations ---"
  grep -iE 'undefined|LaTeX Warning: Reference|Citation.*undefined' /tmp/${f}_c3.log | head -6
done
echo "===== confirm featured figures embedded ====="
grep -oE '(rhc_forest_both|sim_negtransfer_C1_p100|sim_negtransfer_C1_p50)\.pdf' main.log supplemental.log 2>/dev/null | sort | uniq -c

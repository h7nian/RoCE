#!/bin/bash
#SBATCH --job-name=compile_tex
#SBATCH --time=00:20:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=2
#SBATCH --nice=0
#SBATCH --output=diagnosis/face_probe/validation/face_compile_tex_%j.out
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
    echo "[FAIL] $f.tex did not produce a PDF"
  fi
  echo "--- errors/warnings (last pass) ---"
  grep -iE '^!|! LaTeX Error|Undefined control|Runaway|Missing|Overfull \\\\hbox \(.*too wide' /tmp/${f}_c3.log | head -8
  echo "--- undefined refs/citations ---"
  grep -iE 'undefined|LaTeX Warning: Reference|Citation.*undefined' /tmp/${f}_c3.log | head -6
done

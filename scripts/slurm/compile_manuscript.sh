#!/bin/bash
#SBATCH --job-name=roce_tex
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_tex.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_tex.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
DOCS_ROOT="${PROJECT_ROOT}/docs"

module load texlive/2025
cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

for mirrored_file in \
  main.tex \
  supplemental.tex \
  roce_common.sty \
  ref.bib; do
  cmp "${DOCS_ROOT}/${mirrored_file}" \
    "${PROJECT_ROOT}/overleaf/${mirrored_file}"
done

cd "${DOCS_ROOT}"
for document in main supplemental; do
  pdflatex -interaction=nonstopmode -halt-on-error "${document}.tex"
  bibtex "${document}"
  pdflatex -interaction=nonstopmode -halt-on-error "${document}.tex"
  pdflatex -interaction=nonstopmode -halt-on-error "${document}.tex"
  test -s "${document}.pdf"
done

if rg -n \
    'undefined references|Citation .* undefined|Reference .* undefined|There were undefined references' \
    main.log supplemental.log; then
  echo "Unresolved manuscript reference or citation detected." >&2
  exit 1
fi

echo "Compiled main.tex and supplemental.tex without unresolved references."

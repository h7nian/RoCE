#!/bin/bash
#SBATCH --job-name=tex_audit
#SBATCH --time=00:20:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=2
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/face_tex_audit_%j.out
module load texlive/2025 2>/dev/null
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"/docs
echo "=== compile main.tex ==="; pdflatex -interaction=nonstopmode -halt-on-error main.tex >/tmp/main_compile.log 2>&1 && echo "main.tex OK ($(pdfinfo main.pdf 2>/dev/null | grep Pages))" || { echo "main.tex FAILED"; grep -iE '^!|error' /tmp/main_compile.log | head; }
echo "=== compile supplemental.tex ==="; pdflatex -interaction=nonstopmode -halt-on-error supplemental.tex >/tmp/supp_compile.log 2>&1 && echo "supplemental.tex OK" || { echo "supplemental FAILED"; grep -iE '^!|error' /tmp/supp_compile.log | head; }
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
echo "=== audit + key aggregation tests ==="
Rscript -e 'library(testthat); library(RoCE); for(f in c("test-method-reference-alignment-audit.R")) {cat("---",f,"---\n"); tryCatch(test_file(file.path("tests/testthat",f)), error=function(e)cat("ERR:",conditionMessage(e),"\n"))}' 2>&1 | tail -20

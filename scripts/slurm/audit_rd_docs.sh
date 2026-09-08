#!/bin/bash
#SBATCH --job-name=roce_rd_audit
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rd_audit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rd_audit.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p results/direct_tate_mc500_b5000/logs

Rscript -e '
  rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
  if (length(rd_files) == 0L) {
    stop("no generated Rd files found under man/", call. = FALSE)
  }
  options(warn = 1)
  parse_warnings <- character(0)
  invisible(lapply(rd_files, function(path) {
    withCallingHandlers(
      tools::parse_Rd(path),
      warning = function(condition) {
        parse_warnings <<- c(
          parse_warnings,
          sprintf("%s: %s", path, conditionMessage(condition))
        )
        invokeRestart("muffleWarning")
      }
    )
  }))
  if (length(parse_warnings) > 0L) {
    stop(
      sprintf(
        "%d Rd parser warning(s):\n%s",
        length(parse_warnings), paste(unique(parse_warnings), collapse = "\n")
      ),
      call. = FALSE
    )
  }
  message(sprintf("[pass] parsed %d Rd files without warnings", length(rd_files)))
'

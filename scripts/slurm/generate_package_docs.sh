#!/bin/bash
#SBATCH --job-name=roce_roxygen
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_roxygen.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_roxygen.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p results/direct_tate_mc500_b5000/logs

Rscript -e '
  if (!requireNamespace("roxygen2", quietly = TRUE)) {
    stop("roxygen2 is required to regenerate package documentation.",
         call. = FALSE)
  }
  roxygen2::roxygenise(".", roclets = "rd")
  rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
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
  message(sprintf(
    "[pass] regenerated and parsed %d Rd files without warnings",
    length(rd_files)
  ))
'

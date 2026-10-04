#!/usr/bin/env Rscript
# Two separate worker waves fit the same p100 training task. Memory-only
# caches disappear with a worker; the opt-in shared cache survives the wave.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: check_shared_cache_cost.R LIBRARY INPUT_RDS NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
stopifnot(.Platform$OS.type != "windows")
output <- args[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
folds <- readRDS(args[2L])$folds
invisible(RoCE:::set_nuisance_solver_cpp("proximal_newton"))
messages <- RoCE:::.prepare_source_calibration_messages(folds$target_folds,
  folds$source_folds$s1, "s1", 3:10, "one_round", 1L, "binomial", 5, 100L, "min",
  calibration_recipe = "score_derivative")
writeLines(c("1000/site, p100, 200 features, ten original folds, 100 lambdas",
  "Fixed source1/arm1 inner calibration on folds 3:10, final tolerance 1e-10",
  "Separate worker waves; this isolates reuse across outer-fold worker lifetimes",
  "Not a full TATE or coverage experiment"), file.path(output, "scope.txt"))
reference <- NULL
timings <- list()
for (backend in c("memory", "shared")) {
  cache <- RoCE:::.new_nuisance_cache(TRUE, if (backend == "shared") output else NULL)
  for (wave in 1:2) {
    started <- proc.time()[["elapsed"]]
    worker <- parallel::mcparallel({
      fitted <- RoCE:::.fit_source_calibration(folds$source_folds$s1, "s1", 3:10,
        function(site, fold) messages[[paste0("k2_", fold)]], 1L, 1L, 1L, 5,
        100L, 10000L, "min", fit_cache = cache, layout = "compact", tol = 1e-10,
        calibration_recipe = "score_derivative")
      list(fit = fitted, pid = Sys.getpid())
    }, mc.set.seed = FALSE)
    result <- parallel::mccollect(worker)[[1L]]
    if (inherits(result, "try-error")) stop(as.character(result))
    name <- paste(backend, wave, sep = "_")
    saveRDS(result, file.path(output, paste0(name, ".rds")))
    fit <- result$fit
    if (is.null(reference)) reference <- fit
    for (field in c("weight", "outcome")) {
      stopifnot(identical(as.numeric(fit[[field]]), as.numeric(reference[[field]])),
                identical(attr(fit[[field]], "lambda_used"), attr(reference[[field]], "lambda_used")))
    }
    timings[[name]] <- data.frame(backend = backend, wave = wave, pid = result$pid,
      seconds = proc.time()[["elapsed"]] - started,
      cv_seconds = sum(fit$timing[grepl("_cv$", names(fit$timing))]))
    write.csv(do.call(rbind, timings), file.path(output, "timing.csv"), row.names = FALSE)
    cat(name, "completed in", timings[[name]]$seconds, "seconds\n")
  }
  RoCE:::.release_nuisance_cache(cache)
}
stopifnot(length(unique(vapply(timings, function(row) row$pid, integer(1L)))) == 4L,
          timings$shared_2$cv_seconds == 0)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("PASSED", file.path(output, "status.txt"))

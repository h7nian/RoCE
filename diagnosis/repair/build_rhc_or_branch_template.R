#!/usr/bin/env Rscript
# Add a milder shift satisfying the radius3 population moment-feasibility
# prerequisite. The original severe law and all its results stay unchanged.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
original_path <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
source("diagnosis/repair/rhc_model_validation.R")
template <- readRDS(original_path)
original_populations <- template$populations
population <- template$populations$O_case_mix
target <- population$sites$t
row_hashes <- function(x) vapply(seq_len(nrow(x)), function(i)
  digest::digest(as.numeric(x[i, ]), algo = "sha256"), character(1L))
profiles <- row_hashes(template$features)
stopifnot(!anyDuplicated(profiles))
design <- cbind(1, template$features)
decomposition <- qr(design)
stopifnot(decomposition$rank == ncol(design))
diagnostics <- list()
for (site_name in setdiff(names(population$sites), "t")) {
  original_site <- template$observed_structure[[site_name]]
  index <- match(row_hashes(original_site$W_outcome), profiles)
  stopifnot(!anyNA(index))
  own_probability <- tabulate(index, nrow(design)) / length(index)
  site <- population$sites[[site_name]]
  site$probability <- .10 * own_probability + .90 * target$probability
  # Keep the recorded bounded logistic propensity and common true outcome.
  # Only the source feature distribution changes.
  site$outcome <- target$outcome
  population$sites[[site_name]] <- site
  joint <- site$probability * cbind(1 - site$propensity, site$propensity)
  log_weight <- log(target$probability / joint)
  linear_residual <- qr.resid(decomposition, -log_weight)
  for (arm in 0:1) {
    q <- exp(log_weight[, arm + 1L])
    stopifnot(all(q > exp(-3)), all(q < exp(3)),
      max(abs(linear_residual[, arm + 1L])) > 1e-3)
    for (h in list(rep(1, nrow(design)), target$outcome[, arm + 1L] * (1 - target$outcome[, arm + 1L])))
      stopifnot(max(abs(colSums(design * (joint[, arm + 1L] * q * h)) -
        colSums(design * (target$probability * h)))) < 1e-12)
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(site = site_name, arm = arm,
      minimum_true_weight = min(q), maximum_true_weight = max(q),
      radius3_log_margin = 3 - max(abs(log_weight[, arm + 1L])),
      max_log_weight_linear_residual = max(abs(linear_residual[, arm + 1L])))
  }
}
population$description <- paste("Correct common logistic OR; each source feature law is",
  "90% target law plus10% original source empirical law; joint weight is not log-linear",
  "andtrue weights lie strictly inside radius3 bounds")
stopifnot(identical(population$sites$t, original_populations$O_case_mix$sites$t))
template$populations$O_moderate_mix <- population
for (name in names(original_populations)) stopifnot(identical(template$populations[[name]], original_populations[[name]]))
rhc_check_population(population, design)
dir.create(output, recursive = TRUE)
saveRDS(template, file.path(output, "template.rds"))
truth <- read.csv(file.path(dirname(original_path), "population_truth.csv"), stringsAsFactors = FALSE)
new_truth <- truth[truth$scenario == "O_case_mix", ]; new_truth$scenario <- "O_moderate_mix"
write.csv(rbind(truth, new_truth), file.path(output, "population_truth.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(output, "population_model_checks.csv"), row.names = FALSE)
paths <- c(original_path, file.path(output, "template.rds"), "diagnosis/repair/rhc_model_validation.R",
           "diagnosis/repair/build_rhc_or_branch_template.R", "diagnosis/repair/rhc_validation_helpers.R")
jsonlite::write_json(setNames(lapply(paths, function(path) digest::digest(file=path, algo="sha256")),
  normalizePath(paths)), file.path(output, "checked_files.json"), auto_unbox=TRUE, pretty=TRUE)
writeLines(c("The target law,outcome functions andpropensities are unchanged from O_case_mix.",
  "True q is strictly inside radius3 bounds andmatches all initial/h-weighted moments; the shared design has full column rank.",
  "Nonzero log-weight projection residual proves the working log-linear joint-weight model is misspecified.",
  "This verifies key model andmoment-feasibility prerequisites,not every nuisance-rate assumption."),
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, diagnostics))

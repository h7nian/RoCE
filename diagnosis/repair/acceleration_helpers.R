# Shared inspection for numerical acceleration checks; no model fitting.
collect_calibration_models <- function(fit) {
  models <- list()
  add_calibration <- function(calibrated, prefix) {
    stopifnot(is.list(calibrated), is.numeric(calibrated$outcome),
              is.numeric(calibrated$weight), length(calibrated$initial_outcome) > 0L,
              length(calibrated$initial_weight) > 0L)
    models[[paste0(prefix, "/outcome")]] <<- calibrated$outcome
    models[[paste0(prefix, "/weight")]] <<- calibrated$weight
    for (type in c("initial_outcome", "initial_weight")) for (fold in names(calibrated[[type]])) {
      models[[paste(prefix, type, fold, sep = "/")]] <<- calibrated[[type]][[fold]]
    }
  }
  stopifnot(identical(fit$crossfit_levels, 3L),
            setequal(names(fit$arm_results), c("mu1", "mu0")))
  for (arm in names(fit$arm_results)) for (outer in seq_len(fit$n_folds)) {
    block <- fit$arm_results[[arm]]$fold_results[[outer]]
    prefix <- paste(arm, outer, sep = "/")
    add_calibration(block$target_only$calibration, paste0(prefix, "/target/outer"))
    for (inner in names(block$target_only_inner)) {
      add_calibration(block$target_only_inner[[inner]]$calibration, paste(prefix, "target", inner, sep = "/"))
    }
    for (site in names(block$source_results)) {
      source <- block$source_results[[site]]
      models[[paste(prefix, site, "outer/outcome", sep = "/")]] <- source$alpha_ts
      models[[paste(prefix, site, "outer/weight", sep = "/")]] <- source$gamma_s
      for (type in c("per_k2_alpha", "per_k2_gamma")) for (inner in names(source[[type]])) {
        models[[paste(prefix, site, type, inner, sep = "/")]] <- source[[type]][[inner]]
      }
      for (inner in names(source$inner_calibrated)) {
        add_calibration(source$inner_calibrated[[inner]], paste(prefix, site, inner, sep = "/"))
      }
    }
  }
  models
}

compare_calibration_models <- function(reference, candidate, lambda_tolerance = 1e-10) {
  reference <- collect_calibration_models(reference)
  candidate <- collect_calibration_models(candidate)
  stopifnot(length(reference) > 0L, identical(names(reference), names(candidate)))
  do.call(rbind, lapply(names(reference), function(key) {
    a <- candidate[[key]]
    b <- reference[[key]]
    lambda_a <- attr(a, "lambda_used")
    lambda_b <- attr(b, "lambda_used")
    stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)),
              length(lambda_a) == 1L, length(lambda_b) == 1L,
              identical(attr(a, "training_folds"), attr(b, "training_folds")),
              identical(attr(a, "cv_seed"), attr(b, "cv_seed")))
    lambda_equal <- (is.na(lambda_a) && is.na(lambda_b)) ||
      (is.finite(lambda_a) && is.finite(lambda_b) &&
       abs(lambda_a - lambda_b) <= lambda_tolerance * max(1, abs(lambda_b)))
    data.frame(model = key, reference_lambda = lambda_b, candidate_lambda = lambda_a,
      lambda_equal = lambda_equal,
      max_coefficient_difference = max(abs(as.numeric(a) - as.numeric(b))))
  }))
}

#!/usr/bin/env Rscript
# Baseline comparison coverage probe.
#
# Diagnosis-only. For a given config, runs all comparison-method baselines
# (run_all_comparisons: target_only, sample_size, inverse_variance, federated_dr,
# pooled_dr, tilted_aipw) plus the RoCE aggregated estimator (run_crossfit), and
# records each method's bias, se, and CI coverage vs the superpopulation truth.
# Used to show which baselines have bias/coverage problems across C1-C4.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(sprintf("expected 7 args: tag config n_total K p seed n_folds; got %d", length(args)),
       call. = FALSE)
}
tag <- args[[1L]]; config <- args[[2L]]
n_total <- as.integer(args[[3L]]); K_sites <- as.integer(args[[4L]])
p <- as.integer(args[[5L]]); seed <- as.integer(args[[6L]]); n_folds <- as.integer(args[[7L]])
stopifnot(n_total > 0L, K_sites > 0L, p > 0L, seed > 0L, n_folds >= 3L)
if (!(config %in% c("C1", "C2", "C3", "C4"))) stop("config must be C1/C2/C3/C4.", call. = FALSE)

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) library(RoCE)
  else if (requireNamespace("devtools", quietly = FALSE)) devtools::load_all(".")
  else stop("Neither installed RoCE nor devtools available.", call. = FALSE)
})
roce_function <- function(name) {
  ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) stop(sprintf("missing '%s'.", name), call. = FALSE)
  get(name, envir = ns, inherits = FALSE)
}
for (fn in c("generate_simulation_data", "split_data_by_site", "build_crossfit_folds",
             "run_crossfit", "run_all_comparisons")) {
  if (!exists(fn, mode = "function")) assign(fn, roce_function(fn))
}

out_root <- Sys.getenv("C2_CMP_OUTPUT_ROOT", file.path("diagnosis", "c2", "comparison_coverage"))
out_dir <- file.path(out_root, "results"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf("%s_%s_n%d_K%d_p%d_seed%d_kf%d",
                                     tag, config, n_total, K_sites, p, seed, n_folds))

set.seed(seed)
data <- generate_simulation_data(
  n_total = n_total, K = K_sites, p = p, config = config,
  estimand_type = "superpopulation", site_allocation = "model",
  transform_type = "mild", outcome_type = "binary",
  heterogeneity_type = "none", shift_strength = 0.5,
  dgp_type = "roce", warn_ignored = FALSE
)
split <- split_data_by_site(data)
folds <- build_crossfit_folds(split, n_folds)
truth <- as.numeric(data$mu1_true)

get_se <- function(res) {
  if (!is.null(res$se) && is.finite(as.numeric(res$se))) return(as.numeric(res$se))
  v <- suppressWarnings(as.numeric(res$variance))
  if (length(v) && is.finite(v) && v >= 0) return(sqrt(v))
  NA_real_
}
row_for <- function(method, est, se) {
  est <- as.numeric(est); se <- as.numeric(se)
  covered <- if (is.finite(est) && is.finite(se)) abs(est - truth) <= 1.96 * se else NA
  data.frame(tag = tag, config = config, n_total = n_total, K = K_sites, p = p, seed = seed,
             method = method, estimate = est, truth = truth, bias = est - truth,
             se = se, covered = covered, stringsAsFactors = FALSE)
}

rows <- list()
# RoCE aggregated estimator (the proposed method)
set.seed(seed + 200000L)
roce_ok <- tryCatch({
  res <- run_crossfit(split, n_folds = n_folds, communication_mode = "one_round",
                      lambda_selection = "cv", lambda_rule = "min", verbose = FALSE,
                      n_cores = 1L, family = "binomial", A_val = 1L,
                      use_lambda_cache = TRUE, precomputed_folds = folds,
                      nuisance_lambda_rule = "min")
  rows[[length(rows) + 1L]] <<- row_for("face_c", res$estimate, res$se); TRUE
}, error = function(e) { rows[[length(rows) + 1L]] <<- row_for("face_c", NA, NA); FALSE })

# Baselines
set.seed(seed + 400000L)
cmp <- tryCatch(
  run_all_comparisons(split, use_crossfit = TRUE, n_folds = n_folds, family = "binomial", A_val = 1L),
  error = function(e) { cat("[cmp] run_all_comparisons error:", conditionMessage(e), "\n"); NULL })
if (!is.null(cmp)) {
  for (nm in names(cmp)) {
    r <- cmp[[nm]]
    rows[[length(rows) + 1L]] <- row_for(nm, r$estimate, get_se(r))
  }
}

df <- do.call(rbind, rows)
write.csv(df, paste0(prefix, "_comparison.csv"), row.names = FALSE)
cat(sprintf("[cmp-cov] %s config=%s n=%d K=%d roce_ok=%s wrote %d method rows\n",
            tag, config, n_total, K_sites, roce_ok, nrow(df)))
print(df[, c("method", "bias", "se", "covered")], row.names = FALSE)
cat("[cmp-cov] DONE\n")

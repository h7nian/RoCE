#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: profile_p100_stages.R STAGE")
}
stage <- args[[1L]]
valid_stages <- c("face_mu1", "naive_tate", "dr_tate")
if (!stage %in% valid_stages) {
  stop("STAGE must be one of: ", paste(valid_stages, collapse = ", "))
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))

output_root <- Sys.getenv(
  "ROCE_PROFILE_ROOT",
  "results/direct_tate_mc500_b5000/profile_p100"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_root, paste0(stage, ".csv"))
if (file.exists(output_path)) {
  message("[skip] profile already exists: ", output_path)
  quit(save = "no", status = 0L)
}

sim_id <- 1L
config <- "C3"
p <- 100L
K <- 4L
n_site <- 1000L
n_folds <- 5L
nlambda_init <- 100L
allocated_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
source_cores <- min(K, allocated_cores)
if (is.na(source_cores) || source_cores < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.")
}

message(sprintf(
  "[%s] stage=%s sim=%d C=%s p=%d K=%d n_site=%d folds=%d cores=%d",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage, sim_id, config,
  p, K, n_site, n_folds, source_cores
))

set.seed(sim_id)
generated <- generate_simulation_data(
  n_total = n_site * (K + 1L),
  K = K,
  p = p,
  config = config,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  dgp_type = "face",
  ate_deviation = 0,
  n_deviated_sites = 0L,
  warn_ignored = FALSE
)
data_split <- split_data_by_site(generated)
truth_mu1 <- generated$mu1_true
truth_tate <- generated$mu1_true - generated$mu0_true

stage_started <- proc.time()[["elapsed"]]
if (stage == "face_mu1") {
  # build_crossfit_folds() is intentionally package-internal; qualify it so
  # this installed-package profiling script does not depend on load_all().
  folds <- RoCE:::build_crossfit_folds(data_split, n_folds)
  fit <- run_crossfit(
    data_split = data_split,
    n_folds = n_folds,
    communication_mode = "one_round",
    lambda_selection = RoCE:::AGG_WALD_LAMBDA,
    verbose = TRUE,
    n_cores = source_cores,
    nlambda_init = nlambda_init,
    family = "binomial",
    A_val = 1L,
    precomputed_folds = folds
  )
  rows <- data.frame(
    method = "one_round_crossfit_mu1",
    estimate = fit$estimate,
    se = fit$se,
    truth = truth_mu1,
    stringsAsFactors = FALSE
  )
} else {
  methods <- if (stage == "naive_tate") {
    c("sample_size", "inverse_variance")
  } else {
    c("federated_dr", "pooled_dr")
  }
  message("[stage] fitting methods: ", paste(methods, collapse = ", "))
  fits <- run_all_comparisons_tate(
    data_split = data_split,
    family = "binomial",
    n_folds = n_folds,
    variance_method = "bootstrap",
    include_tilted = FALSE,
    methods = methods,
    n_cores = source_cores
  )
  rows <- do.call(rbind, lapply(methods, function(method) {
    data.frame(
      method = paste0(method, "_ate"),
      estimate = fits[[method]]$estimate,
      se = fits[[method]]$se,
      truth = truth_tate,
      stringsAsFactors = FALSE
    )
  }))
}

elapsed <- proc.time()[["elapsed"]] - stage_started
rows$bias <- rows$estimate - rows$truth
rows$stage <- stage
rows$elapsed_seconds <- elapsed
rows$sim_id <- sim_id
rows$config <- config
rows$p <- p
rows$K <- K
rows$rho <- 0
rows$n_site <- n_site
rows$n_folds <- n_folds
rows$nlambda_init <- nlambda_init
rows$allocated_cores <- allocated_cores

temporary_path <- tempfile(
  pattern = paste0(".", stage, "_"),
  tmpdir = output_root,
  fileext = ".csv"
)
write.csv(rows, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path)
  stop("failed to atomically move profile output to ", output_path)
}
message(sprintf("[done] %s elapsed %.1f seconds; wrote %s",
                stage, elapsed, output_path))

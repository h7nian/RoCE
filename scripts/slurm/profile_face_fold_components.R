#!/usr/bin/env Rscript

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))

output_root <- Sys.getenv(
  "ROCE_PROFILE_ROOT",
  "results/direct_tate_mc500_b5000/profile_face_fold_components"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

read_integer_env <- function(name, default, lower, upper = Inf) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower || value > upper) {
    stop(sprintf("%s must be an integer in [%s, %s].", name, lower, upper))
  }
  value
}

read_numeric_env <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    stop(sprintf("%s must be finite and lie in [%s, %s].", name, lower, upper))
  }
  value
}

sim_id <- read_integer_env("ROCE_PROFILE_SIM_ID", 1L, 1L)
config <- Sys.getenv("ROCE_PROFILE_CONFIG", "C3")
if (!config %in% c("C1", "C2", "C3")) {
  stop("ROCE_PROFILE_CONFIG must be one of C1, C2, or C3.")
}
p <- 100L
K <- read_integer_env("ROCE_PROFILE_K", 1L, 1L)
n_site <- read_integer_env("ROCE_PROFILE_N_SITE", 1000L, 10L)
n_folds <- 5L
k1 <- read_integer_env("ROCE_PROFILE_OUTER_FOLD", 1L, 1L, n_folds)
A_val <- read_integer_env("ROCE_PROFILE_ARM", 1L, 0L, 1L)
source_index <- read_integer_env("ROCE_PROFILE_SOURCE_INDEX", 1L, 1L, K)
rho <- read_numeric_env("ROCE_PROFILE_RHO", 0, lower = 0)
nuisance_nlambda <- read_integer_env(
  "ROCE_NUISANCE_NLAMBDA", 100L, 2L
)
target_nlambda <- read_integer_env(
  "ROCE_TARGET_NLAMBDA", nuisance_nlambda, 2L
)
nuisance_max_iter <- read_integer_env(
  "ROCE_NUISANCE_MAX_ITER", RoCE:::MAX_ITER_DEFAULT, 1L
)

format_filename_number <- function(value) {
  gsub("\\.", "p", format(value, scientific = FALSE, trim = TRUE))
}
arm_label <- if (A_val == 1L) "mu1" else "mu0"
output_path <- file.path(
  output_root,
  sprintf(
    paste0(
      "%s_sim%d_p%d_K%d_rho%s_fold%d_source%d_%s_",
      "targetgrid%d_sourcegrid%d_iter%d.csv"
    ),
    config, sim_id, p, K, format_filename_number(rho), k1, source_index,
    arm_label, target_nlambda, nuisance_nlambda, nuisance_max_iter
  )
)
if (file.exists(output_path)) {
  message("[skip] profile already exists: ", output_path)
  quit(save = "no", status = 0L)
}

set.seed(sim_id)
generated <- generate_simulation_data(
  n_total = n_site * (K + 1L),
  K = K,
  p = p,
  config = config,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  dgp_type = "face",
  ate_deviation = rho,
  n_deviated_sites = if (rho > 0) 1L else 0L,
  warn_ignored = FALSE
)
data_split <- split_data_by_site(generated)
folds <- RoCE:::build_crossfit_folds(data_split, n_folds)
target_folds <- folds$target_folds
source_folds <- folds$source_folds
source_site <- names(source_folds)[[source_index]]
secondary_folds <- setdiff(seq_len(n_folds), k1)

target_state <- RoCE:::with_seed(
  RoCE:::.crossfit_work_seed("one_round", A_val),
  {
    target_initial_models <- list()
    target_summaries <- list()
    cached_lambda <- NULL
    target_started_at <- proc.time()[["elapsed"]]
    for (k2 in secondary_folds) {
      key <- paste0("k2_", k2)
      training_folds <- setdiff(seq_len(n_folds), c(k1, k2))
      target_train <- RoCE:::combine_folds(target_folds, training_folds)
      target_calibration <- RoCE:::materialize_fold(target_folds, k2)
      alpha_initial <- fit_initial_outcome(
        target_train$W_outcome,
        target_train$Y,
        target_train$A,
        A_val = A_val,
        lambda = cached_lambda,
        nlambda = target_nlambda,
        family = "binomial"
      )
      if (is.null(cached_lambda)) {
        selected_lambda <- attr(alpha_initial, "lambda_used")
        if (!is.null(selected_lambda) && is.finite(selected_lambda)) {
          cached_lambda <- selected_lambda
        }
      }
      target_initial_models[[key]] <- alpha_initial
      target_summaries[[key]] <- list(
        mean_phi = c(1, colMeans(target_train$Z_site)),
        mean_grad_psi_init = RoCE:::.mean_glm_gradient_site_basis(
          target_calibration$W_outcome,
          target_calibration$Z_site,
          alpha_initial,
          RoCE:::FAMILY_BINOMIAL,
          RoCE:::LINK_LOGIT
        )
      )
    }
    list(
      initial_models = target_initial_models,
      summaries = target_summaries,
      elapsed_seconds = proc.time()[["elapsed"]] - target_started_at
    )
  }
)
target_initial_models <- target_state$initial_models
target_summaries <- target_state$summaries
target_initial_seconds <- target_state$elapsed_seconds

get_fold_inputs <- function(site, k2) {
  key <- paste0("k2_", k2)
  list(
    mean_phi = target_summaries[[key]]$mean_phi,
    mean_grad_psi_init = target_summaries[[key]]$mean_grad_psi_init,
    alpha_init = target_initial_models[[key]]
  )
}

source_result <- RoCE:::with_seed(
  RoCE:::.crossfit_work_seed(
    "one_round", A_val, k1, source_index
  ),
  RoCE:::process_source_site(
    s = source_site,
    source_folds = source_folds,
    target_folds = target_folds,
    k1 = k1,
    n_folds = n_folds,
    A_val = A_val,
    M_tau = RoCE:::M_TAU_DEFAULT,
    M_tau_inference = RoCE:::M_TAU_INFERENCE_DEFAULT,
    data_split = data_split,
    combine_cache = new.env(hash = TRUE, parent = emptyenv()),
    get_fold_inputs = get_fold_inputs,
    family_int = RoCE:::FAMILY_BINOMIAL,
    link_int = RoCE:::LINK_LOGIT,
    use_lambda_cache = TRUE,
    nuisance_nlambda = nuisance_nlambda,
    nuisance_max_iter = nuisance_max_iter
  )
)

row <- as.data.frame(as.list(c(
  target_initial_seconds = as.numeric(target_initial_seconds),
  source_result$timing,
  RoCE:::.summarize_source_nuisance_fit_diagnostics(
    setNames(list(source_result), source_site)
  )
)))
row$sim_id <- sim_id
row$config <- config
row$p <- p
row$K <- K
row$source_index <- source_index
row$source_site <- source_site
row$rho <- rho
row$n_site <- n_site
row$n_folds <- n_folds
row$target_nlambda <- target_nlambda
row$nuisance_nlambda <- nuisance_nlambda
row$nuisance_max_iter <- nuisance_max_iter
row$outer_fold <- k1
row$A_val <- A_val
row$estimate <- source_result$mu_ts
row$estimate_finite <- is.finite(source_result$mu_ts)
row$initial_density_ratio_lambda <- mean(
  source_result$nuisance_fit_diagnostics$initial_density_ratio_lambdas
)
row$calibrated_density_ratio_lambda <-
  source_result$nuisance_fit_diagnostics$calibrated_density_ratio_lambda
row$calibrated_density_ratio_update_ratio <-
  source_result$nuisance_fit_diagnostics$calibrated_density_ratio_update_ratio
row$calibrated_density_ratio_max_abs_coefficient <-
  source_result$nuisance_fit_diagnostics$
    calibrated_density_ratio_max_abs_coefficient
row$calibrated_outcome_lambda <-
  source_result$nuisance_fit_diagnostics$calibrated_outcome_lambda

temporary_path <- tempfile(
  pattern = ".face_fold_components_",
  tmpdir = output_root,
  fileext = ".csv"
)
write.csv(row, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path)
  stop("failed to atomically move profile output to ", output_path)
}
message("[done] wrote ", output_path)

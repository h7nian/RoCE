# Published RHC profile. Use run.R rather than sourcing this file alone.
collation <- Sys.setlocale("LC_COLLATE", "C.UTF-8")
if (!identical(collation, "C.UTF-8")) stop("RHC reproduction requires the C.UTF-8 locale (Linux/WSL).")
list(
  seed = 42L, n_bootstrap = 5000L,
  baselines = c("sample_size", "inverse_variance", "federated_dr", "pooled_dr"),
  build_data = function() {
    recode <- stats::setNames(rep(NA_character_, 2L), c("Medicare & Medicaid", "No insurance"))
    cohort_args <- list(outcome = "death30", site_var = "ninsclas", site_recode = recode,
                        covariate_profile = "historical")
    split_args <- list(K = 4L, target_site = "Private", phi = base::identity, seed = 42L)
    cohort <- do.call(RoCE::build_rhc_cohort, cohort_args)
    data <- do.call(RoCE::build_rhc_data_split, c(list(cohort = cohort), split_args))
    stopifnot(nrow(cohort) == 5039L, ncol(data$t$W_outcome) == 61L)
    folds <- RoCE:::build_crossfit_folds(data, 10L)
    prepared <- RoCE:::.prepare_rhc_outer_preprocessing(cohort, data, folds, cohort_args, split_args)
    list(data_split = data, folds = prepared$folds, comparison_data = prepared$baseline_data,
      metadata = list(target_site = "Private", site_mapping = attr(data, "site_mapping"),
        preprocessing = "outer_fold", inner_preprocessing = "Shared outer-training transformation",
        communication = "Federated emulation with shared covariate preprocessing"))
  },
  fit = list(n_folds = 10L, communication_mode = "one_round", family = "binomial",
    crossfit_layers = 2L, aggregation_mode = "joint_tate", lambda_selection = .5,
    target_nuisance_method = "hou_calibrated", source_validation_method = "outer_fit",
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "logistic",
      source_nuisance_method = "calibrated", target_radius = 5),
    calibration_layout = "compact", nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    nuisance_cv_certificate = TRUE, M_tau = 3, M_tau_inference = 3,
    nlambda_init = 100L, nuisance_lambda_rule = "min", verbose = TRUE)
)

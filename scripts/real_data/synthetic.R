# Small end-to-end installation check; deliberately NOT a reportable analysis.
list(
  seed = 20261003L, n_bootstrap = 200L,
  baselines = c("sample_size", "inverse_variance", "federated_dr", "pooled_dr"),
  build_data = function() {
    generated <- RoCE::generate_simulation_data(n_total = 3000L, K = 2L, p = 4L,
      config = "C1", dgp_type = "bounded", outcome_type = "binary")
    data <- RoCE::split_data_by_site(generated)
    for (site in names(data)) {
      colnames(data[[site]]$W_outcome) <- colnames(data[[site]]$Z_site) <- paste0("x", 1:4)
    }
    list(data_split = data, comparison_data = data,
      folds = RoCE:::build_crossfit_folds(data, 4L), metadata = list(purpose = "Synthetic installation check"))
  },
  fit = list(n_folds = 4L, communication_mode = "one_round", family = "binomial",
    crossfit_layers = 2L, aggregation_mode = "joint_tate", lambda_selection = .5,
    target_nuisance_method = "hou_calibrated", source_validation_method = "outer_fit",
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "logistic",
      source_nuisance_method = "calibrated", target_radius = 5),
    calibration_layout = "compact", nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    nuisance_cv_certificate = TRUE, M_tau = 3, M_tau_inference = 3,
    nlambda_init = 5L, nuisance_lambda_rule = "min", verbose = FALSE)
)

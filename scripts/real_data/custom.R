# Copy this template and edit the scientific specification for your data.
# ROCE_DATA_RDS points to a local data.frame, e.g. inside an All of Us workspace.
list(
  seed = 42L, n_bootstrap = 5000L,
  baselines = c("sample_size", "inverse_variance", "federated_dr", "pooled_dr"),
  build_data = function() {
    path <- Sys.getenv("ROCE_DATA_RDS")
    if (!nzchar(path)) stop("Set ROCE_DATA_RDS to your local data-frame RDS file.")
    frame <- readRDS(path)
    data <- site_data_from_frame(frame, site_column = "site", treatment_column = "A",
      outcome_column = "Y", feature_columns = c("x1", "x2", "x3", "x4"), target_site = "target")
    folds <- RoCE:::build_crossfit_folds(data, 10L)
    list(data_split = data, folds = folds, comparison_data = data,
      metadata = list(site_mapping = attr(data, "site_mapping"),
        preprocessing = "Pre-specified numeric features; no sample-fitted transformations"))
  },
  # Starting configuration, NOT settings validated for a different population.
  fit = list(n_folds = 10L, communication_mode = "one_round", family = "binomial",
    crossfit_layers = 2L, aggregation_mode = "joint_tate", lambda_selection = .5,
    target_nuisance_method = "hou_calibrated", source_validation_method = "outer_fit",
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "logistic",
      source_nuisance_method = "calibrated", target_radius = 5),
    calibration_layout = "compact", nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    nuisance_cv_certificate = TRUE, M_tau = 3, M_tau_inference = 3,
    nlambda_init = 100L, nuisance_lambda_rule = "min", verbose = TRUE)
)

# Optional RHC preprocessing learned without an outer evaluation fold.
# Internal calibration/CV reuses this outer-training transformation. This mode
# isolates the outer-preprocessing boundary; it is not nested CV preprocessing.

.validate_rhc_preprocessing_rows <- function(rows, n) {
  if (is.null(rows)) return(invisible(NULL))
  if (!is.numeric(rows) || !length(rows) || anyNA(rows) ||
      any(!is.finite(rows)) || any(rows != as.integer(rows)) ||
      any(rows < 1L | rows > n) || anyDuplicated(rows)) {
    stop("preprocessing_rows must identify distinct retained-cohort training rows.", call. = FALSE)
  }
  invisible(NULL)
}

# Fixed labels in the historical C.UTF-8 order, not estimated from a fold.
.rhc_categorical_levels <- function(profile = "historical") {
  result <- list(sex = c("Female", "Male"), race = c("black", "other", "white"),
    income = c("> $50k", "$11-$25k", "$25-$50k", "Under $11k"),
    ninsclas = c("Medicaid", "Medicare", "Medicare & Medicaid", "No insurance",
      "Private", "Private & Medicare"), ca = c("Metastatic", "No", "Yes"),
    dnr1 = c("No", "Yes"), cat1 = c("ARF", "CHF", "Cirrhosis", "Colon Cancer",
      "Coma", "COPD", "Lung Cancer", "MOSF w/Malignancy", "MOSF w/Sepsis"))
  if (profile == "grouped_log_missing") {
    result$cat1 <- c(setdiff(result$cat1, c("Colon Cancer", "Lung Cancer")), "Solid Cancer")
  }
  result
}

.rhc_training_preprocessed_data <- function(training_rows, cohort_arguments, split_arguments) {
  cohort_arguments$preprocessing_rows <- training_rows
  cohort_arguments$categorical_levels <- .rhc_categorical_levels(cohort_arguments$covariate_profile)
  cohort <- do.call(build_rhc_cohort, cohort_arguments)
  split_arguments$cohort <- cohort
  split_arguments$preprocessing_rows <- training_rows
  # build_rhc_data_split retains its historical RNG initialization. Isolate it
  # here so preprocessing cannot change nuisance CV or bootstrap draws.
  data <- with_seed(split_arguments$seed, do.call(build_rhc_data_split, split_arguments))
  attr(data, "training_preprocessing") <- list(training_rows = as.integer(training_rows),
    feature_center = attr(data, "feature_center"), feature_scale = attr(data, "feature_scale"),
    continuous_fill = vapply(attr(cohort, "continuous_vars"), function(variable)
      stats::median(cohort[[variable]][training_rows]), numeric(1L)),
    categorical_levels = cohort_arguments$categorical_levels)
  data
}

.rhc_cohort_rows_by_site <- function(cohort, data) {
  components <- attr(data, "site_components")
  rows <- lapply(names(data), function(site) {
    key <- if (site == "t") "target" else site
    index <- which(cohort$site_var %in% components[[key]])
    if (length(index) != data[[site]]$n ||
        !identical(cohort$A[index], data[[site]]$A) ||
        !identical(cohort$Y[index], data[[site]]$Y)) {
      stop("RHC preprocessing must preserve the cohort and site row order.", call. = FALSE)
    }
    index
  })
  names(rows) <- names(data)
  if (!identical(sort(unlist(rows, use.names = FALSE)), seq_len(nrow(cohort)))) {
    stop("RHC site rows must partition the retained cohort.", call. = FALSE)
  }
  rows
}

.rhc_replace_fold_data <- function(folds, data) {
  result <- lapply(folds, function(fold) fold)
  attr(result, ".data_ref") <- data
  result
}

.rhc_baseline_preprocessor <- function(site, site_rows, n_total, cohort_arguments, split_arguments) {
  force(site); force(site_rows); force(n_total)
  force(cohort_arguments); force(split_arguments)
  function(training_rows) {
    .validate_rhc_preprocessing_rows(training_rows, length(site_rows))
    excluded <- site_rows[setdiff(seq_along(site_rows), training_rows)]
    data <- .rhc_training_preprocessed_data(setdiff(seq_len(n_total), excluded),
      cohort_arguments, split_arguments)
    data[[site]][c("W_outcome", "Z_site")]
  }
}

.prepare_rhc_outer_preprocessing <- function(cohort, data, folds, cohort_arguments, split_arguments) {
  if (!identical(split_arguments$phi, base::identity)) {
    stop("Outer-fold RHC preprocessing currently requires phi=identity.", call. = FALSE)
  }
  if (is.null(cohort_arguments$raw)) cohort_arguments$raw <- load_rhc_raw()
  rows <- .rhc_cohort_rows_by_site(cohort, data)
  site_folds <- c(list(t = folds$target_folds), folds$source_folds)
  n_folds <- length(folds$target_folds)
  prepared <- lapply(names(data), function(site) vector("list", n_folds))
  names(prepared) <- names(data)
  parameters <- vector("list", n_folds)
  for (outer in seq_len(n_folds)) {
    excluded <- unlist(lapply(names(data), function(site)
      rows[[site]][site_folds[[site]][[outer]]$original_idx]), use.names = FALSE)
    training_rows <- setdiff(seq_len(nrow(cohort)), excluded)
    transformed <- .rhc_training_preprocessed_data(training_rows, cohort_arguments, split_arguments)
    for (site in names(data)) {
      if (!identical(colnames(transformed[[site]]$W_outcome), colnames(data[[site]]$W_outcome)) ||
          !identical(transformed[[site]]$A, data[[site]]$A) ||
          !identical(transformed[[site]]$Y, data[[site]]$Y)) {
        stop("Outer preprocessing changed the fixed RHC feature dictionary or records.", call. = FALSE)
      }
      prepared[[site]][[outer]] <- .rhc_replace_fold_data(site_folds[[site]], transformed[[site]])
    }
    parameters[[outer]] <- attr(transformed, "training_preprocessing")
  }
  attr(folds$target_folds, ".outer_preprocessed_views") <- prepared$t
  for (site in names(folds$source_folds)) {
    attr(folds$source_folds[[site]], ".outer_preprocessed_views") <- prepared[[site]]
  }
  attr(folds, "preprocessing") <- list(mode = "outer_fold", parameters = parameters,
    inner_scope = "Initial and calibration CV reuse the outer-training transformation",
    communication_scope = "Pooled covariate preprocessing within outer training; federated emulation")
  baseline_data <- data
  for (site in names(data)) {
    attr(baseline_data[[site]], "training_preprocessor") <- .rhc_baseline_preprocessor(
      site, rows[[site]], nrow(cohort), cohort_arguments, split_arguments)
  }
  list(folds = folds, baseline_data = baseline_data)
}

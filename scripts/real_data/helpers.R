# Numeric, pre-specified features only: no whole-cohort learned preprocessing.
site_data_from_frame <- function(data, site_column, treatment_column, outcome_column,
                                 feature_columns, target_site) {
  roles <- c(site_column, treatment_column, outcome_column)
  if (!is.data.frame(data) || anyDuplicated(names(data)) || length(roles) != 3L ||
      anyDuplicated(roles) || !length(feature_columns) || anyDuplicated(feature_columns) ||
      any(feature_columns %in% roles) || !all(c(roles, feature_columns) %in% names(data))) {
    stop("Supply distinct site/treatment/outcome columns and existing, distinct feature columns.")
  }
  labels <- as.character(data[[site_column]])
  if (anyNA(labels) || any(!nzchar(labels)) || length(target_site) != 1L ||
      is.na(target_site) || !target_site %in% labels) stop("Invalid site labels or target_site.")
  features <- data[feature_columns]
  if (!all(vapply(features, is.numeric, logical(1L)))) {
    stop("Features must be numeric with a common pre-specified encoding across sites.")
  }
  if (!is.numeric(data[[treatment_column]]) || !is.numeric(data[[outcome_column]])) {
    stop("Treatment and outcome must be numeric; treatment must be coded 0/1.")
  }
  feature_matrix <- as.matrix(features)
  if (any(!is.finite(feature_matrix))) stop("Features must be finite; no automatic imputation is performed.")
  source_labels <- unique(labels[labels != target_site])
  if (!length(source_labels)) stop("At least one source site is required.")
  ordered_labels <- c(target_site, source_labels)
  result <- lapply(ordered_labels, function(label) {
    rows <- which(labels == label)
    x <- feature_matrix[rows, , drop = FALSE]
    list(X = x, X_dagger = x, W_outcome = x, Z_site = x,
         A = data[[treatment_column]][rows], Y = data[[outcome_column]][rows], n = length(rows))
  })
  names(result) <- c("t", paste0("s", seq_along(source_labels)))
  attr(result, "site_mapping") <- stats::setNames(ordered_labels, names(result))
  result
}

real_data_method_row <- function(label, fit) {
  stopifnot(length(fit$estimate) == 1L, length(fit$se) == 1L,
            is.finite(fit$estimate), is.finite(fit$se), fit$se > 0)
  data.frame(method = label, estimate = fit$estimate, se = fit$se,
    ci_lower = fit$estimate - stats::qnorm(.975) * fit$se,
    ci_upper = fit$estimate + stats::qnorm(.975) * fit$se)
}

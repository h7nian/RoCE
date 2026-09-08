# First-order sampling-design decomposition for equal-within-arm fold allocation.
# Input is the observation-mass gradient of the complete layer under study.
# This is not a scalar SE correction and does not differentiate nuisance fits.
.balanced_arm_fold_gradient <- function(gradient, treatment, fold) {
  n <- length(gradient)
  if (!is.numeric(gradient) || !is.null(dim(gradient)) || n < 2L ||
      any(!is.finite(gradient)) || length(treatment) != n || length(fold) != n ||
      anyNA(treatment) || any(!treatment %in% c(0, 1)) ||
      !is.numeric(fold) || any(!is.finite(fold)) || any(fold != floor(fold)) || any(fold < 1)) {
    stop("invalid gradient, treatment or fold records")
  }
  labels <- sort(unique(fold))
  if (!identical(as.integer(labels), seq_along(labels))) stop("fold labels must be consecutive")
  cell_counts <- table(factor(fold, levels = labels), factor(treatment, levels = 0:1))
  if (any(cell_counts == 0) || any(apply(cell_counts, 2L, function(x) diff(range(x))) > 1)) {
    stop("requires nonempty, equally allocated treatment-arm folds, up to integer remainders")
  }
  within <- gradient
  cell_means <- matrix(NA_real_, length(labels), 2L)
  for (k in labels) for (arm in 0:1) {
    selected <- fold == k & treatment == arm
    value <- mean(gradient[selected])
    cell_means[k, arm+1L] <- value
    within[selected] <- gradient[selected]-value
  }
  # A site treatment-total fluctuation is allocated equally to every fold,
  # not as independent fold-specific treatment-proportion fluctuations.
  treatment_slope <- n*mean(cell_means[, 2L]-cell_means[, 1L])
  probability <- mean(treatment)
  composition <- treatment_slope*(treatment-probability)/n
  adjusted <- within+composition
  stopifnot(abs(sum(adjusted)) < 1e-10*max(1, sum(abs(adjusted))),
            abs(sum(within*composition)) < 1e-12*max(1, sum(adjusted^2)))
  list(gradient = adjusted, within_cell = within, treatment_composition = composition,
       treatment_slope = treatment_slope, cell_counts = cell_counts,
       cell_means = cell_means, variance = sum(adjusted^2),
       integer_remainders_differentiated = FALSE)
}

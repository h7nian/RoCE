#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: analyze_full_target_tails.R FULL_FITS_RDS OUTPUT_CSV")
artifact <- readRDS(args[1L])
target <- artifact$data_split$t
variants <- names(artifact$variants$mu1[[1L]]$outer)
rows <- list()
for (variant in variants) for (A_val in c(1L, 0L)) {
  arm_name <- if (A_val == 1L) "mu1" else "mu0"
  probability <- prediction <- pseudo_outcome <- rep(NA_real_, target$n)
  assignment_count <- integer(target$n)
  for (fold in seq_along(artifact$folds$target_folds)) {
    indices <- artifact$folds$target_folds[[fold]]$original_idx
    fit <- artifact$variants[[arm_name]][[fold]]$outer[[variant]]
    probability[indices] <- if (A_val == 1L) fit$prop_scores else 1 - fit$prop_scores
    prediction[indices] <- fit$m_pred
    pseudo_outcome[indices] <- fit$estimate + fit$varphi_ot
    assignment_count[indices] <- assignment_count[indices] + 1L
  }
  stopifnot(all(assignment_count == 1L), all(is.finite(pseudo_outcome)))
  selected <- which(target$A == A_val)
  worst <- selected[which.max(1 / probability[selected])]
  true_probability <- if (A_val == 1L) artifact$truth$propensity[worst] else
    1 - artifact$truth$propensity[worst]
  rows[[length(rows) + 1L]] <- data.frame(
    variant, A_val, max_observed_inverse_probability = max(1 / probability[selected]),
    n_observed_probability_below_005 = sum(probability[selected] < .05),
    worst_row = worst, worst_fitted_probability = probability[worst],
    worst_true_probability = true_probability, worst_outcome = target$Y[worst],
    worst_outcome_prediction = prediction[worst],
    worst_pseudo_outcome = pseudo_outcome[worst], estimate = mean(pseudo_outcome)
  )
}
rows <- do.call(rbind, rows)
write.csv(rows, args[2L], row.names = FALSE)
print(rows, row.names = FALSE)

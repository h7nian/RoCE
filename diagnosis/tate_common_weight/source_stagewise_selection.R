# Requires target_projection_validation.R for the shared complete-fold selector.
.select_source_stagewise_candidates <- function(losses) {
  candidates <- c("null", "c0.5", "c1", "c2")
  required <- c("final_candidate", "initial_candidate", "arm", "fold", "n_validation",
                "final_risk", "initial_risk", "final_failed", "initial_failed")
  if (!is.data.frame(losses) || !all(required %in% names(losses)) || nrow(losses) != 128L ||
      !setequal(losses$final_candidate, candidates) || !setequal(losses$initial_candidate, candidates) ||
      !setequal(losses$arm, 0:1) || !setequal(losses$fold, 2:5) ||
      anyDuplicated(paste(losses$final_candidate, losses$initial_candidate, losses$arm, losses$fold)) ||
      !is.logical(losses$final_failed) || !is.logical(losses$initial_failed) ||
      anyNA(losses$final_failed) || anyNA(losses$initial_failed) ||
      any(!is.finite(losses$n_validation)) || any(losses$n_validation < 1) ||
      any(losses$n_validation != floor(losses$n_validation))) stop("invalid complete source candidate table")
  if (any(!is.finite(losses$final_risk[!losses$final_failed])) ||
      any(!is.finite(losses$initial_risk[!losses$final_failed & !losses$initial_failed]))) stop("nonfinite successful risk")
  for (fold in 2:5) {
    if (length(unique(losses$n_validation[losses$fold == fold])) != 1L) stop("inconsistent validation sample counts")
  }
  for (candidate in candidates) for (arm in 0:1) for (fold in 2:5) {
    x <- losses[losses$final_candidate == candidate & losses$arm == arm & losses$fold == fold, ]
    if (length(unique(x$final_failed)) != 1L || length(unique(x$final_risk)) != 1L)
      stop("final fit or risk changed across initial branches")
  }
  null_final <- losses$final_candidate == "null"
  if (any(losses$final_failed[null_final]) || any(losses$final_risk[null_final] != 0))
    stop("null final candidate must succeed with zero risk")
  null_initial <- losses$initial_candidate == "null" & !losses$final_failed
  if (any(losses$initial_failed[null_initial]) || any(losses$initial_risk[null_initial] != 0))
    stop("null initial candidate must succeed with zero risk")
  collapse_arms <- function(x, candidate_column, risk_column, failure) {
    rows <- list()
    for (candidate in candidates) for (fold in 2:5) {
      selected <- x[[candidate_column]] == candidate & x$fold == fold
      pair <- x[selected, ]
      stopifnot(nrow(pair) == 2L, setequal(pair$arm, 0:1))
      failed <- any(failure[selected])
      rows[[length(rows)+1L]] <- data.frame(candidate, fold, n_validation = pair$n_validation[1],
        loss = if (failed) NA_real_ else sum(pair[[risk_column]]), failed)
    }
    do.call(rbind, rows)
  }
  final_table <- losses[losses$initial_candidate == "null", ]
  final_losses <- collapse_arms(final_table, "final_candidate", "final_risk", final_table$final_failed)
  final_selection <- .select_target_projection_candidate(final_losses)
  branch <- losses[losses$final_candidate == final_selection$selected, ]
  initial_losses <- collapse_arms(branch, "initial_candidate", "initial_risk", branch$final_failed | branch$initial_failed)
  initial_selection <- .select_target_projection_candidate(initial_losses)
  if (final_selection$selected == "null") {
    # With zero downstream coefficients every initial RHS vanishes.
    if (any(initial_losses$failed) || any(initial_losses$loss != 0)) stop("null final branch must have zero initial risks")
    stopifnot(initial_selection$selected == "null")
  }
  list(final_selection = final_selection, initial_selection = initial_selection,
       final_losses = final_losses, initial_losses = initial_losses)
}

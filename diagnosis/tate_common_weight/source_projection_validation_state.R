# Validate a fitted source equation on new data while keeping its fold weights.
.source_projection_validation_state <- function(state, target, source) {
  stopifnot(length(state$keys) > 0L, nrow(target$W) > 0L, nrow(source$W) > 0L)
  for (site in c("target", "source")) {
    heldout <- if (site == "target") target else source
    stopifnot(!any(heldout$index %in% state[[site]]$index))
    for (piece in state$pieces) for (part in c("training", "calibration"))
      stopifnot(!any(heldout$index %in% piece[[paste(site, part, sep = "_")]]$index))
  }
  fractions <- vapply(state$pieces, function(piece)
    nrow(piece$source_calibration$W)/nrow(state$source$W), numeric(1L))
  stopifnot(identical(names(fractions), state$keys), all(fractions > 0), abs(sum(fractions)-1) < 1e-12)
  validation <- state
  validation$target <- target
  validation$source <- source
  validation$pieces <- setNames(lapply(state$keys, function(key)
    list(target_training = target, source_training = source,
         target_calibration = target, source_calibration = source)), state$keys)
  # Each repeated validation block contributes fraction * mean(moment),
  # not an additional full mean. Target summaries retain equal fold weights.
  validation$source_calibration_denominators <- as.list(nrow(source$W)/fractions)
  validation
}

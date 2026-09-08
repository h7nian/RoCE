roce_parse_semicolon_numeric_grid <- function(
    value, name, require_zero_first = FALSE, require_positive = FALSE) {
  if (length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop(name, " must be one nonempty semicolon-separated grid.",
         call. = FALSE)
  }
  parsed <- suppressWarnings(as.numeric(
    strsplit(as.character(value), ";", fixed = TRUE)[[1L]]
  ))
  invalid <- length(parsed) < 1L || any(!is.finite(parsed)) ||
    anyDuplicated(parsed) || any(parsed < 0) ||
    (require_positive && any(parsed <= 0))
  if (invalid) {
    qualifier <- if (require_positive) "positive" else "nonnegative"
    stop(
      name, " must contain unique finite ", qualifier,
      " values separated by semicolons.", call. = FALSE
    )
  }
  if (require_zero_first &&
      (length(parsed) < 2L || !identical(parsed[[1L]], 0))) {
    stop(name, " must start at zero and contain at least two values.",
         call. = FALSE)
  }
  parsed
}

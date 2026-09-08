# Runtime recording for the isolated, sequential nuisance-refit diagnostic.
# No installed package files are changed. The traced function is restored on
# success and failure, and the recorder never samples random numbers.

roce_cv_partition_record_valid <- function(record) {
  if (!is.list(record) ||
      !all(c("caller", "n_folds", "cv_group_id", "cv_fold_id") %in% names(record))) return(FALSE)
  groups <- record$cv_group_id
  n_folds <- record$n_folds
  if (!is.character(record$caller) || length(record$caller) != 1L ||
      is.na(record$caller) || !nzchar(record$caller) ||
      !typeof(groups) %in% c("integer", "double") ||
      !is.null(dim(groups)) || is.object(groups) || !length(groups) ||
      any(!is.finite(groups)) || any(groups < 1) || any(groups != floor(groups)) ||
      any(groups > .Machine$integer.max) ||
      !is.numeric(n_folds) || length(n_folds) != 1L || !is.finite(n_folds) ||
      n_folds < 2 || n_folds != floor(n_folds) || n_folds > length(groups)) return(FALSE)
  folds <- record$cv_fold_id
  if (is.null(folds)) return(anyDuplicated(groups) == 0L)
  if (!typeof(folds) %in% c("integer", "double") ||
      !is.null(dim(folds)) || is.object(folds) || length(folds) != length(groups) ||
      any(!is.finite(folds)) || any(folds != floor(folds)) ||
      any(folds < 1) || any(folds > n_folds)) return(FALSE)
  counts <- tabulate(as.integer(folds), nbins = as.integer(n_folds))
  if (any(counts == 0L) || any(counts == length(groups))) return(FALSE)
  all(vapply(split(folds, groups), function(x) length(unique(x)) == 1L, logical(1L)))
}

roce_with_nuisance_cv_audit <- function(fit_function) {
  if (!is.function(fit_function)) stop("CV audit requires a zero-argument fit function.")
  namespace <- asNamespace("RoCE")
  symbol <- ".make_nuisance_cv_fold_id"
  if (!exists(symbol, envir = namespace, inherits = FALSE)) {
    stop("Origin-grouped refitting requires the tested grouped-CV RoCE package.")
  }
  if (inherits(get(symbol, envir = namespace), "functionWithTrace")) {
    stop("Refusing to replace an existing nuisance-CV trace.")
  }
  state <- new.env(parent = emptyenv())
  state$records <- list()
  record_partition <- function(groups, n_folds, caller, folds) {
    record <- list(caller = caller, n_folds = n_folds,
                   cv_group_id = groups, cv_fold_id = folds)
    # Do not mask a fitting error if the function exits before building folds.
    record$valid <- tryCatch(roce_cv_partition_record_valid(record),
                             error = function(e) FALSE)
    state$records[[length(state$records) + 1L]] <- record
    invisible(NULL)
  }
  exit_expression <- substitute(
    recorder(cv_group_id, n_folds, caller, returnValue()),
    list(recorder = record_partition)
  )
  trace(symbol, exit = exit_expression, print = FALSE, where = namespace)
  on.exit(untrace(symbol, where = namespace), add = TRUE)
  fitted <- tryCatch(fit_function(), error = function(e) e)
  list(fitted = fitted, partitions = state$records)
}

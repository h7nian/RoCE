# Optional grouped cross-validation support for nuisance-model tuning.

.validate_nuisance_cv_group_id <- function(cv_group_id, n, caller) {
  if (is.null(cv_group_id)) return(NULL)
  if (!is.character(caller) || length(caller) != 1L || is.na(caller) ||
      !nzchar(caller)) {
    stop("caller must be one nonempty string.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop(caller, ": n must be one positive integer.", call. = FALSE)
  }
  if (!typeof(cv_group_id) %in% c("integer", "double") || !is.null(dim(cv_group_id)) ||
      is.object(cv_group_id) || length(cv_group_id) != n ||
      anyNA(cv_group_id) || any(!is.finite(cv_group_id)) ||
      any(cv_group_id < 1) || any(cv_group_id != floor(cv_group_id)) ||
      any(cv_group_id > .Machine$integer.max)) {
    stop(
      caller,
      ": cv_group_id must be a positive finite integer vector of length ",
      as.integer(n), ".",
      call. = FALSE
    )
  }
  as.integer(cv_group_id)
}

.make_nuisance_cv_fold_id <- function(cv_group_id, n_folds, caller) {
  if (is.null(cv_group_id)) return(NULL)
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, length(cv_group_id), caller
  )
  if (!is.numeric(n_folds) || length(n_folds) != 1L || is.na(n_folds) ||
      !is.finite(n_folds) || n_folds < 2 || n_folds != floor(n_folds) ||
      n_folds > .Machine$integer.max) {
    stop(caller, ": n_folds must be one integer >= 2.", call. = FALSE)
  }
  n_folds <- as.integer(n_folds)
  group_levels <- unique(cv_group_id)
  # Preserve the exact legacy glmnet path when every row is its own group.
  # In particular, do not call sample() or otherwise advance the RNG.
  if (length(group_levels) == length(cv_group_id)) return(NULL)
  if (length(group_levels) < n_folds) {
    stop(
      caller, ": grouped nuisance CV has only ", length(group_levels),
      " unique groups for ", n_folds, " folds.", call. = FALSE
    )
  }
  group_fold <- sample(rep(
    seq_len(n_folds), length.out = length(group_levels)
  ))
  as.integer(group_fold[match(cv_group_id, group_levels)])
}

# C++ selectors take explicit fold IDs for arm-filtered rows. Keep NULL/unique
# groups on their original internal random-fold path, without consuming RNG.
.call_nuisance_cv_with_groups <- function(cv_function, args, cv_group_id,
                                          A, A_val, n_folds, caller) {
  cv_fold_id <- .make_nuisance_cv_fold_id(
    if (is.null(cv_group_id)) NULL else cv_group_id[A == A_val],
    n_folds, caller
  )
  if (!is.null(cv_fold_id)) args$cv_fold_id <- cv_fold_id
  do.call(cv_function, args)
}

# Group metadata must not be silently lost through stale precomputed views.
# Repeated origins also must not cross outer folds: callers supplying a
# bootstrap multiset must retain its original group-consistent partition.
.validate_nuisance_cv_fold_views <- function(site_data, folds, caller) {
  groups <- .validate_nuisance_cv_group_id(site_data$cv_group_id, site_data$n, caller)
  view_groups <- .validate_nuisance_cv_group_id(
    attr(folds, ".data_ref")$cv_group_id, site_data$n, caller
  )
  if (!identical(groups, view_groups)) {
    stop(caller, ": precomputed fold views have different nuisance CV group IDs.", call. = FALSE)
  }
  if (is.null(groups)) return(invisible(TRUE))
  if (!is.list(folds) || length(folds) < 2L) stop(caller, ": invalid grouped fold views.")
  indices <- lapply(folds, `[[`, "original_idx")
  if (any(vapply(indices, function(idx) {
    !is.numeric(idx) || !length(idx) || any(!is.finite(idx)) ||
      any(idx != floor(idx)) || any(idx < 1) || any(idx > site_data$n)
  }, logical(1L))) ||
      !identical(sort(as.integer(unlist(indices))), seq_len(site_data$n))) {
    stop(caller, ": grouped fold views must partition all row positions once.", call. = FALSE)
  }
  memberships <- split(rep(seq_along(folds), lengths(indices)),
                        groups[unlist(indices, use.names = FALSE)])
  if (any(vapply(memberships, function(x) length(unique(x)) > 1L, logical(1L)))) {
    stop(caller, ": an origin crosses outer folds; supply group-consistent precomputed_folds.", call. = FALSE)
  }
  invisible(TRUE)
}

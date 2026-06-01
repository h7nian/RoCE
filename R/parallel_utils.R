# parallel_utils.R - Parallelization utilities for FACE-HD cross-fitting
#
# Provides setup_parallel() and parallel_lapply() used by cross_fitting_algorithms.R
# to process multiple source sites in parallel.

# =============================================================================
# PARALLELIZATION UTILITIES
# =============================================================================

#' Setup parallel backend
#'
#' @param n_cores Number of cores to use. If NULL, uses sequential processing.
#'   If -1, uses all available cores minus 1.
#' @return Number of cores to use (1 for sequential)
setup_parallel <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    return(1)
  }
  if (length(n_cores) != 1L || !is.numeric(n_cores) || is.na(n_cores) ||
      abs(n_cores - round(n_cores)) > sqrt(.Machine$double.eps)) {
    stop("setup_parallel: n_cores must be NULL, -1, or a positive integer.",
         call. = FALSE)
  }
  if (n_cores == 1) {
    return(1)
  }
  if (n_cores == 0 || n_cores < -1) {
    stop("setup_parallel: n_cores must be NULL, -1, or a positive integer.",
         call. = FALSE)
  }
  
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop(sprintf(
      "setup_parallel: requested n_cores=%s but the 'parallel' namespace is not available.",
      format(n_cores)
    ), call. = FALSE)
  }
  
  max_cores <- parallel::detectCores()
  if (length(max_cores) != 1L || is.na(max_cores) || max_cores < 1L) {
    stop("setup_parallel: parallel::detectCores() did not return a positive core count.",
         call. = FALSE)
  }
  
  if (n_cores == -1) {
    # Use all cores minus 1
    n_cores <- max(1, max_cores - 1)
  } else if (n_cores > max_cores) {
    warning(sprintf(
      "setup_parallel: requested n_cores=%d exceeds detected cores=%d; using %d cores.",
      as.integer(n_cores), as.integer(max_cores), as.integer(max_cores)
    ), call. = FALSE)
    n_cores <- max_cores
  }
  
  return(as.integer(max(1, n_cores)))
}

#' Apply function over list with optional parallelization
#'
#' @param x List to iterate over
#' @param fun Function to apply
#' @param n_cores Number of cores (1 for sequential)
#' @param ... Additional arguments to fun
#' @return List of results
parallel_lapply <- function(x, fun, n_cores = 1, ...) {
  # Skip parallelism when there is at most one element: fork/cluster setup
  # overhead exceeds the compute savings for trivially small inputs.  This
  # also covers the K=1 source-site case where no work can be distributed.
  if (n_cores == 1 || length(x) <= 1) {
    return(lapply(x, fun, ...))
  } else if (.Platform$OS.type == "windows") {
    # mclapply doesn't work on Windows; use parLapply via socket cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Export only the functions transitively needed by 'fun' to workers.
    # Data variables are captured via fun's closure automatically.
    if (!requireNamespace("codetools", quietly = TRUE)) {
      stop("parallel_lapply: codetools is required for Windows socket workers; install codetools or run with n_cores = 1.",
           call. = FALSE)
    }
    needed_fns <- tryCatch({
      available_fns <- Filter(function(nm) is.function(get(nm, envir = globalenv())),
                              ls(envir = globalenv()))
      # Recursively find all global functions referenced by fun and its callees
      seen <- character(0)
      queue <- intersect(codetools::findGlobals(fun, merge = FALSE)$functions, available_fns)
      while (length(queue) > 0) {
        fn_name <- queue[1L]
        queue <- queue[-1L]
        if (fn_name %in% seen) next
        seen <- c(seen, fn_name)
        sub_refs <- tryCatch(
          codetools::findGlobals(get(fn_name, envir = globalenv()), merge = FALSE)$functions,
          error = function(e) {
            stop(sprintf("parallel_lapply: codetools::findGlobals failed for '%s': %s",
                         fn_name, conditionMessage(e)), call. = FALSE)
          }
        )
        queue <- c(queue, setdiff(intersect(sub_refs, available_fns), seen))
      }
      seen
    }, error = function(e) {
      stop(sprintf("parallel_lapply: static dependency analysis failed; refusing to export all global functions automatically. Original error: %s",
                   conditionMessage(e)), call. = FALSE)
    })
    if (length(needed_fns) > 0) {
      parallel::clusterExport(cl, needed_fns, envir = globalenv())
    }
    return(parallel::parLapply(cl, x, fun, ...))
  } else {
    # Use mclapply for Unix-like systems (fork-based, more efficient)
    return(parallel::mclapply(x, fun, ..., mc.cores = n_cores))
  }
}

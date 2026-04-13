# parallel_utils.R - Parallelization utilities for FACE-C cross-fitting
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
  if (is.null(n_cores) || n_cores == 1) {
    return(1)
  }
  
  # Load parallel package
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("parallel package not available, using sequential processing")
    return(1)
  }
  
  max_cores <- parallel::detectCores()
  
  if (n_cores == -1) {
    # Use all cores minus 1
    n_cores <- max(1, max_cores - 1)
  } else {
    n_cores <- min(n_cores, max_cores)
  }
  
  return(max(1, n_cores))
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
    needed_fns <- tryCatch({
      if (!requireNamespace("codetools", quietly = TRUE)) {
        stop("codetools not available")
      }
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
            warning(sprintf("parallel_lapply: codetools::findGlobals failed for '%s': %s. Skipping transitive deps.",
                            fn_name, conditionMessage(e)))
            character(0)
          }
        )
        queue <- c(queue, setdiff(intersect(sub_refs, available_fns), seen))
      }
      seen
    }, error = function(e) {
      # Fallback: export all global functions if static analysis fails
      warning(sprintf("parallel_lapply: static dependency analysis failed: %s. Exporting all global functions.",
                      conditionMessage(e)))
      Filter(function(nm) is.function(get(nm, envir = globalenv())),
             ls(envir = globalenv()))
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

# Initial nuisance fits are functions of their own training subset and tuning
# inputs. A per-run cache can reuse an identical fit; it cannot share a selected
# lambda or a warm start between different training problems.

NUISANCE_TRAINING_POLICY <- "training_subset_v3_calibration"

.validate_calibration_control <- function(control = NULL, target_nuisance_method = "hou_calibrated") {
  defaults <- list(recipe = "legacy", target_propensity_initialization = "logistic",
                   target_radius = NULL, source_nuisance_method = "calibrated")
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control)) || any(!names(control) %in% c(names(defaults), "source_lambda_rules"))))) {
    stop(paste0("calibration_control must be a named list with recipe, ",
      "target_propensity_initialization, target_radius, source_nuisance_method and/or source_lambda_rules."), call. = FALSE)
  }
  for (name in names(control)) defaults[name] <- control[name]
  for (name in c("recipe", "target_propensity_initialization", "source_nuisance_method")) {
    if (!is.character(defaults[[name]]) || length(defaults[[name]]) != 1L || is.na(defaults[[name]])) {
      stop("calibration_control$", name, " must be one character value.", call. = FALSE)
    }
  }
  defaults$recipe <- match.arg(defaults$recipe, c("legacy", "score_derivative"))
  defaults$target_propensity_initialization <- match.arg(
    defaults$target_propensity_initialization, c("logistic", "calibrated"))
  defaults$source_nuisance_method <- match.arg(defaults$source_nuisance_method, c("calibrated", "standard"))
  if (!is.null(defaults$source_lambda_rules)) {
    .resolve_source_lambda_rules(defaults$source_lambda_rules, "min")
    if (defaults$recipe != "score_derivative" || defaults$source_nuisance_method != "calibrated") {
      stop("source_lambda_rules requires score_derivative calibrated sources.", call. = FALSE)
    }
  } else {
    defaults$source_lambda_rules <- NULL
  }
  if (!is.null(defaults$target_radius)) {
    radius <- defaults$target_radius
    if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || radius <= 0) {
      stop("calibration_control$target_radius must be positive or NULL.", call. = FALSE)
    }
    defaults$target_radius <- as.numeric(radius)
  }
  if (target_nuisance_method == "lasso" &&
      (defaults$target_propensity_initialization != "logistic" || !is.null(defaults$target_radius))) {
    stop("Target calibration settings require target_nuisance_method='hou_calibrated'.", call. = FALSE)
  }
  defaults
}

# Override individual source fitting stages while target fits and initial OR
# messages continue to use nuisance_lambda_rule. Missing stages inherit it.
.resolve_source_lambda_rules <- function(rules = NULL, lambda_rule = "min") {
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, ".resolve_source_lambda_rules")
  result <- list(initial_weight = lambda_rule, weight = lambda_rule, outcome = lambda_rule)
  if (is.null(rules)) return(result)
  if (!is.list(rules) || !length(rules) || is.null(names(rules)) || anyNA(names(rules)) ||
      anyDuplicated(names(rules)) || any(!names(rules) %in% names(result))) {
    stop("source_lambda_rules must name initial_weight, weight and/or outcome.", call. = FALSE)
  }
  for (stage in names(rules)) {
    value <- rules[[stage]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% c("min", "1se")) {
      stop("source_lambda_rules$", stage, " must be min or 1se.", call. = FALSE)
    }
    result[[stage]] <- value
  }
  result
}

.validate_score_calibration_radius <- function(radius) {
  # Below this bound the declared symmetric tilt truncation precedes every
  # additional numerical/ratio guard in the source evaluation score.
  maximum <- min(log(WEIGHT_MAX), -log(WEIGHT_MIN), log(INFERENCE_RATIO_MAX))
  if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) ||
      radius <= 0 || radius > maximum) {
    stop("Score-derivative calibration requires a finite positive fitting radius <= ",
         signif(maximum, 8), ".", call. = FALSE)
  }
  invisible(radius)
}

.validate_logical_control <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.set_nuisance_solver <- function(nuisance_solver) {
  if (is.null(nuisance_solver)) return(NULL)
  solver <- match.arg(nuisance_solver, c("proximal_newton", "coordinate_descent"))
  set_nuisance_solver_cpp(solver)
}

.restore_nuisance_solver <- function(previous_solver) {
  if (!is.null(previous_solver)) set_nuisance_solver_cpp(previous_solver)
  invisible(NULL)
}

.nuisance_cv_certificate_enabled <- function() {
  .validate_logical_control(getOption("RoCE.nuisance_cv_certificate", FALSE),
                            "nuisance_cv_certificate")
}

.set_nuisance_cv_certificate <- function(value) {
  if (is.null(value)) {
    .nuisance_cv_certificate_enabled()
    return(NULL)
  }
  value <- .validate_logical_control(value, "nuisance_cv_certificate")
  options(RoCE.nuisance_cv_certificate = value)
}

.restore_nuisance_cv_certificate <- function(previous) {
  if (!is.null(previous)) options(previous)
  invisible(NULL)
}

.validate_nuisance_cache_flag <- function(use_lambda_cache, caller) {
  if (!is.logical(use_lambda_cache) || length(use_lambda_cache) != 1L ||
      is.na(use_lambda_cache)) {
    stop(caller, ": use_lambda_cache must be TRUE or FALSE.", call. = FALSE)
  }
  use_lambda_cache
}

.nuisance_training_key <- function(fitter, site, training_folds, A_val) {
  if (!is.character(fitter) || length(fitter) != 1L || is.na(fitter) ||
      !nzchar(fitter) || !is.character(site) || length(site) != 1L ||
      is.na(site) || !nzchar(site) || !is.numeric(training_folds) ||
      !length(training_folds) || anyNA(training_folds) ||
      any(!is.finite(training_folds)) || any(training_folds < 1) ||
      any(training_folds != floor(training_folds)) || anyDuplicated(training_folds)) {
    stop(".nuisance_training_key: invalid fitter, site, or training folds.", call. = FALSE)
  }
  A_val <- .validate_A_val(A_val, ".nuisance_training_key")
  paste0(nchar(fitter, type = "bytes"), ":", fitter, ":",
         nchar(site, type = "bytes"), ":", site, ":", A_val, ":",
         paste(sort(training_folds), collapse = ","))
}

.nuisance_training_seed <- function(key) {
  # Arithmetic stays below 2^53. A seed collision only shares a random stream;
  # cache hits additionally require exact equality of all fitting inputs.
  seed <- 900013
  for (byte in as.integer(charToRaw(enc2utf8(key)))) {
    seed <- (131 * seed + byte) %% 2147483646
  }
  as.integer(seed + 1)
}

.new_nuisance_cache <- function(enabled, directory = NULL, checkpoint_dir = NULL) {
  enabled <- .validate_nuisance_cache_flag(enabled, ".new_nuisance_cache")
  persistent <- !is.null(checkpoint_dir)
  if (persistent && !enabled) stop("checkpoint_dir requires use_lambda_cache=TRUE.", call. = FALSE)
  if (persistent) {
    if (!is.null(directory)) {
      stop("Use only one of nuisance_cache_dir and checkpoint_dir.", call. = FALSE)
    }
    if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
        is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
      stop("checkpoint_dir must be a nonempty directory path.", call. = FALSE)
    }
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    directory <- checkpoint_dir
  }
  if (!is.null(directory) && !enabled) {
    stop("nuisance_cache_dir requires use_lambda_cache=TRUE.", call. = FALSE)
  }
  if (!enabled) return(NULL)
  cache <- new.env(parent = emptyenv())
  if (!is.null(directory)) {
    if (!is.character(directory) || length(directory) != 1L || is.na(directory) ||
        !dir.exists(directory) || file.access(directory, 2L) != 0L) {
      stop("nuisance_cache_dir must name an existing writable directory.", call. = FALSE)
    }
    if (persistent) {
      # Completed fits are immutable, keyed by exact inputs and installed code.
      # A requeued worker replays the same seeds and recovers these fits.
      package_root <- system.file(package = "RoCE")
      code_files <- c(file.path(package_root, "R", "RoCE.rdb"),
        list.files(file.path(package_root, "libs"), pattern = "[.](so|dll)$",
                   recursive = TRUE, full.names = TRUE))
      if (length(code_files) < 2L || any(!file.exists(code_files))) {
        stop("Persistent nuisance checkpoints require an installed RoCE package.", call. = FALSE)
      }
      code_key <- digest::digest(vapply(code_files, function(file)
        digest::digest(file = file, algo = "sha256"), character(1L)), algo = "sha256")
      path <- file.path(normalizePath(directory), paste0("nuisance_", code_key))
      dir.create(path, showWarnings = FALSE)
      if (!dir.exists(path)) stop("Could not create nuisance checkpoint directory.", call. = FALSE)
    } else {
      path <- tempfile("roce_nuisance_", tmpdir = normalizePath(directory))
      if (!dir.create(path)) stop("Could not create the temporary nuisance cache.", call. = FALSE)
      attr(cache, "owner_pid") <- Sys.getpid()
    }
    attr(cache, "shared_directory") <- path
  }
  cache
}

.release_nuisance_cache <- function(cache) {
  path <- attr(cache, "shared_directory")
  if (!is.null(path) && identical(attr(cache, "owner_pid"), Sys.getpid())) {
    unlink(path, recursive = TRUE)
  }
  invisible(NULL)
}

.cached_nuisance_result <- function(fit) {
  if (identical(attr(fit, "lambda_rule"), "degenerate_constant")) .record_or_degenerate_fold()
  for (field in c("cv_seconds", "final_fit_seconds")) {
    if (!is.null(attr(fit, field))) attr(fit, field) <- 0
  }
  fit
}

.fit_nuisance_training_subset <- function(fitter, arguments, site,
                                           training_folds, fit_cache = NULL,
                                           cv_seed_fitter = NULL) {
  key <- .nuisance_training_key(fitter, site, training_folds, arguments$A_val)
  if (!is.null(fit_cache) && !is.environment(fit_cache)) {
    stop(".fit_nuisance_training_subset: fit_cache must be an environment or NULL.",
         call. = FALSE)
  }
  inputs <- list(arguments = arguments, solver = nuisance_solver_cpp(),
                 cv_certificate = .nuisance_cv_certificate_enabled())
  seed_key <- key
  if (!is.null(cv_seed_fitter)) {
    seed_key <- .nuisance_training_key(cv_seed_fitter, site, training_folds, arguments$A_val)
    inputs$cv_seed_fitter <- cv_seed_fitter
  }
  shared_directory <- attr(fit_cache, "shared_directory")
  input_hash <- if (is.null(shared_directory)) NULL else digest::digest(
    list(policy = NUISANCE_TRAINING_POLICY, key = key, inputs = inputs), algo = "sha256")
  entry <- if (is.null(fit_cache)) NULL else fit_cache[[key]]
  if (!is.null(entry)) {
    matches <- if (is.null(shared_directory)) identical(entry$inputs, inputs) else
      identical(entry$input_hash, input_hash)
    if (matches) return(.cached_nuisance_result(entry$fit))
  }
  cache_file <- if (is.null(shared_directory)) NULL else
    file.path(shared_directory, paste0(input_hash, ".rds"))
  if (!is.null(cache_file) && file.exists(cache_file)) {
    entry <- tryCatch(readRDS(cache_file), error = function(error) {
      stop("Could not read shared nuisance cache: ", conditionMessage(error), call. = FALSE)
    })
    if (!is.list(entry) || !identical(entry$input_hash, input_hash) || is.null(entry$fit) ||
        !identical(entry$fit_hash, digest::digest(entry$fit, algo = "sha256"))) {
      stop("Shared nuisance cache failed its input/model integrity check.", call. = FALSE)
    }
    fit_cache[[key]] <- entry
    return(.cached_nuisance_result(entry$fit))
  }
  seed <- .nuisance_training_seed(seed_key)
  fit <- tryCatch(with_seed(seed, do.call(fitter, arguments)), error = function(condition) {
    context <- list(fitter = fitter, site = site, arm = arguments$A_val,
                    training_folds = sort(as.integer(training_folds)), cv_seed = seed)
    stop(structure(list(
      message = sprintf("%s failed at site %s, arm %s, training folds {%s}: %s",
        fitter, site, arguments$A_val, paste(context$training_folds, collapse = ","),
        conditionMessage(condition)), call = NULL, parent = condition, nuisance_context = context
    ), class = c("nuisance_fit_error", "error", "condition")))
  })
  if (!is.null(attr(fit, "converged")) && !isTRUE(attr(fit, "converged"))) {
    stop(sprintf("%s did not converge at site %s on training folds {%s}; lambda=%g, KKT=%g.",
      fitter, site, paste(sort(training_folds), collapse = ","),
      attr(fit, "lambda_used") %||% NA_real_, attr(fit, "kkt_residual") %||% NA_real_), call. = FALSE)
  }
  attr(fit, "cv_seed") <- seed
  attr(fit, "training_folds") <- sort(as.integer(training_folds))
  if (!is.null(fit_cache)) {
    entry <- if (is.null(shared_directory)) list(inputs = inputs, fit = fit) else
      list(input_hash = input_hash, fit = fit, fit_hash = digest::digest(fit, algo = "sha256"))
    fit_cache[[key]] <- entry
    if (!is.null(cache_file)) {
      temporary <- tempfile("fit_", tmpdir = shared_directory)
      on.exit(unlink(temporary), add = TRUE)
      saveRDS(entry, temporary, compress = FALSE)
      if (!file.rename(temporary, cache_file)) {
        stop("Could not atomically store the shared nuisance fit.", call. = FALSE)
      }
    }
  }
  fit
}

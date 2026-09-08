#!/usr/bin/env Rscript

.ritcv_frozen_package <-
  "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7"
.ritcv_frozen_workflow <-
  "9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4"
.ritcv_frozen_manifest <-
  "1dc5b92cf549b417a71eadf1cfe4be50bf8a8219e147e171d53f073ff763a16d"

.ritcv_empty_warning_events <- function() data.frame(
  event_index = integer(), model = character(), arm = integer(),
  k1 = integer(), k2 = integer(), message = character(),
  condition_class = character(), condition_call = character(),
  stringsAsFactors = FALSE
)

.ritcv_output_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || basename(path) %in% c("", ".", "..")) {
    stop("OUTPUT_DIR must be one concrete path.", call. = FALSE)
  }
  output <- file.path(normalizePath(dirname(path), mustWork = FALSE), basename(path))
  if (file.exists(output) || dir.exists(output)) stop("OUTPUT_DIR already exists.", call. = FALSE)
  output
}

.ritcv_integer_scalar <- function(x, name, lower, upper) {
  if (!typeof(x) %in% c("integer", "double") || length(x) != 1L ||
      !is.null(dim(x)) || is.object(x) || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(name, " must be one integer in [", lower, ", ", upper, "].",
         call. = FALSE)
  }
  as.integer(x)
}

.ritcv_build_glmnet_args <- function(
    x, y, family, n_cv_folds, nlambda, maxit,
    nuisance_fold_id = NULL, keep = TRUE) {
  args <- list(x = x, y = y, family = family, alpha = 1,
               nlambda = nlambda, maxit = maxit, keep = keep)
  if (is.null(nuisance_fold_id)) {
    # This is the exact fit_glmnet_cv NULL/all-unique branch. Omitting nfolds
    # would silently select cv.glmnet's default ten folds instead.
    args$nfolds <- n_cv_folds
  } else {
    args$foldid <- nuisance_fold_id
  }
  args
}

.ritcv_prepare_glmnet_args <- function(
    x, y, family, n_cv_folds, nlambda, maxit, cv_group_id, caller,
    fold_maker, keep = TRUE) {
  fold_id <- fold_maker(cv_group_id, n_cv_folds, caller)
  list(
    args = .ritcv_build_glmnet_args(
      x, y, family, n_cv_folds, nlambda, maxit, fold_id, keep
    ),
    fold_id = fold_id
  )
}

.ritcv_validate_prediction <- function(observed, expected, tolerance = 1e-10) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.null(dim(tolerance)) || is.object(tolerance) || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("prediction tolerance must be one positive finite numeric scalar.",
         call. = FALSE)
  }
  valid_vector <- function(x) is.numeric(x) && typeof(x) %in% c("integer", "double") &&
    is.null(dim(x)) && !is.object(x)
  if (!valid_vector(observed) || !valid_vector(expected) ||
      length(observed) != length(expected) || !length(observed) ||
      any(!is.finite(observed)) || any(!is.finite(expected))) {
    stop("replayed and saved predictions must be equal-length finite numeric vectors.",
         call. = FALSE)
  }
  error <- max(abs(observed - expected))
  if (!is.finite(error) || error > tolerance) {
    stop("replayed prediction identity exceeded tolerance: ",
         format(error, digits = 17), call. = FALSE)
  }
  error
}

.ritcv_check_sha <- function(directory, expected_files, sha256_file) {
  path <- file.path(directory, "sha256.txt")
  lines <- readLines(path, warn = FALSE)
  files <- substring(lines, 67L); hashes <- substr(lines, 1L, 64L)
  actual <- list.files(directory, all.files = TRUE, no.. = TRUE)
  if (length(lines) != length(expected_files) ||
      !setequal(files, expected_files) || anyDuplicated(files) ||
      !setequal(actual, c(expected_files, "sha256.txt")) ||
      any(file.info(file.path(directory, actual))$isdir) ||
      any(substr(lines, 65L, 66L) != "  ") ||
      any(!grepl("^[0-9a-f]{64}$", hashes)) ||
      any(vapply(file.path(directory, files), sha256_file,
                 character(1L)) != hashes)) {
    stop("checksum payload validation failed: ", directory, call. = FALSE)
  }
  sha256_file(path)
}

.ritcv_read_gate <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines) || any(!grepl("^[^=]+=[^=]*$", lines))) {
    stop("malformed independent seed audit gate.", call. = FALSE)
  }
  keys <- sub("=.*$", "", lines)
  if (anyDuplicated(keys)) stop("duplicate independent seed audit gate keys.", call. = FALSE)
  stats::setNames(sub("^[^=]*=", "", lines), keys)
}

.ritcv_require_standard_cv <- function(y, model, family) {
  # fit_glmnet_cv uses a constant-outcome fallback before CV in this case.
  # Do not fit an extra CV model and call it a replay of an absent invocation.
  if (identical(model, "OR") && identical(family, "binomial") &&
      min(table(factor(y, levels = c(0, 1)))) < 8L) {
    stop("saved target outcome uses the pre-CV constant fallback; standard CV replay is not applicable.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.ritcv_publish_failure <- function(output, message, state, sha256_file, atomic_writer) {
  atomic_writer(output, function(stage) {
    write.csv(data.frame(
      status = "failed", failure_message = conditionMessage(message),
      replay_completed = FALSE, source_bundle = state$bundle,
      bundle_checksum_manifest_fingerprint = state$bundle_hash,
      seed_audit_checksum_manifest_fingerprint = state$audit_hash,
      task_id = state$task_id, sim_id = state$sim_id,
      current_model = state$current$model, current_arm = state$current$arm,
      current_k1 = state$current$k1, current_k2 = state$current$k2,
      completed_call_count = length(state$diagnostics), stringsAsFactors = FALSE),
              file.path(stage, "attempt_failure.csv"), row.names = FALSE)
    diagnostics <- if (length(state$diagnostics)) do.call(rbind, state$diagnostics) else
      data.frame()
    warnings <- if (length(state$warning_events)) do.call(rbind, state$warning_events) else
      .ritcv_empty_warning_events()
    write.csv(diagnostics, file.path(stage, "completed_calls.csv"), row.names = FALSE)
    write.csv(warnings, file.path(stage, "warning_events.csv"), row.names = FALSE)
    writeLines(c("independent_target_cv_replay=failed", "replay_completed=FALSE"),
               file.path(stage, "metadata.txt"))
    files <- c("attempt_failure.csv", "completed_calls.csv", "warning_events.csv",
               "metadata.txt")
    writeLines(paste(vapply(file.path(stage, files), sha256_file, character(1L)),
                     files, sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "independent target CV replay failure")
}

replay_independent_target_cv_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop("usage: replay_independent_target_cv.R BUNDLE_DIR PACKAGE_LIBRARY OUTPUT_DIR")
  }
  output <- .ritcv_output_path(args[[3L]])
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
  import_files <- c(
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R",
    "diagnosis/tate_common_weight/run_independent_inference_pilot.R",
    "scripts/slurm/result_provenance.R", "scripts/slurm/atomic_output.R"
  )
  script_hash <- roce_sha256_file(
    "diagnosis/tate_common_weight/replay_independent_target_cv.R"
  )
  import_hash <- .weight_calibration_files_fingerprint(
    import_files, roce_sha256_file
  )
  bundle <- normalizePath(args[[1L]], mustWork = TRUE)
  package_library <- normalizePath(args[[2L]], mustWork = TRUE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  state <- new.env(parent = emptyenv())
  state$bundle <- bundle; state$bundle_hash <- NA_character_
  state$audit_hash <- NA_character_
  state$task_id <- NA_integer_; state$sim_id <- NA_integer_
  state$current <- list(model = NA_character_, arm = NA_integer_,
                        k1 = NA_integer_, k2 = NA_integer_)
  state$diagnostics <- list(); state$warning_events <- list()

  attempt <- tryCatch({
    bundle_hash <- .ritcv_check_sha(bundle, c(
      "results.csv", "inference_audit.csv", "diagnostic_qc.csv",
      "artifacts.rds", "metadata.txt"
    ), roce_sha256_file)
    state$bundle_hash <- bundle_hash
    root <- dirname(bundle)
    manifest_path <- file.path(root, "manifest.csv")
    if (roce_sha256_file(manifest_path) != .ritcv_frozen_manifest) {
      stop("ROOT manifest fingerprint is not frozen.")
    }
    manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
    audit <- normalizePath(file.path(root, "audits", basename(bundle)), mustWork = TRUE)
    audit_hash <- .ritcv_check_sha(audit, c(
      "audit_passed.txt", "summary.csv", "numerical_checks.csv",
      "nuisance_diagnostics.csv"
    ), roce_sha256_file)
    state$audit_hash <- audit_hash
    gate <- .ritcv_read_gate(file.path(audit, "audit_passed.txt"))
    required <- c("independent_inference_seed_audit", "bundle_directory",
      "bundle_checksum_manifest_fingerprint", "package_fingerprint",
      "workflow_fingerprint", "manifest_fingerprint", "task_id", "sim_id")
    if (!all(required %in% names(gate)) ||
        gate[["independent_inference_seed_audit"]] != "passed" ||
        normalizePath(gate[["bundle_directory"]], mustWork = TRUE) != bundle ||
        gate[["bundle_checksum_manifest_fingerprint"]] != bundle_hash ||
        gate[["package_fingerprint"]] != .ritcv_frozen_package ||
        gate[["workflow_fingerprint"]] != .ritcv_frozen_workflow ||
        gate[["manifest_fingerprint"]] != .ritcv_frozen_manifest) {
      stop("independent seed audit binding failed.")
    }

    .libPaths(c(package_library, .libPaths()))
    suppressPackageStartupMessages(library(RoCE))
    loaded <- normalizePath(dirname(find.package("RoCE")), mustWork = TRUE)
    actual_package <- .weight_calibration_installed_package_fingerprint(
      package_library, roce_sha256_file
    )
    if (loaded != package_library || actual_package != .ritcv_frozen_package) {
      stop("replay did not load the exact frozen v19 package.")
    }
    saved <- readRDS(file.path(bundle, "artifacts.rds"))
    task_id <- .ritcv_integer_scalar(suppressWarnings(as.numeric(gate[["task_id"]])),
                                     "task_id", 1L, 100L)
    sim_id <- .ritcv_integer_scalar(suppressWarnings(as.numeric(gate[["sim_id"]])),
                                    "sim_id", 10001L, 10100L)
    state$task_id <- task_id; state$sim_id <- sim_id
    manifest_task <- roce_validate_inference_pilot_manifest(manifest, task_id)
    if (manifest_task$sim_id != sim_id || sim_id != 10000L + task_id ||
        saved$task$task_id != task_id || saved$task$sim_id != sim_id) {
      stop("task/simulation identity failed.")
    }
    rows <- read.csv(file.path(bundle, "results.csv"), stringsAsFactors = FALSE)
    required_provenance <- c("task_id", "sim_id", "rho", "method",
      "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint")
    if (!is.data.frame(rows) || nrow(rows) != 84L ||
        !all(required_provenance %in% names(rows)) ||
        !identical(sort(unique(as.numeric(rows$rho))), .inference_pilot_rhos()) ||
        any(rows$task_id != task_id) || any(rows$sim_id != sim_id) ||
        any(rows$package_fingerprint != .ritcv_frozen_package) ||
        any(rows$workflow_fingerprint != .ritcv_frozen_workflow) ||
        any(rows$manifest_fingerprint != .ritcv_frozen_manifest)) {
      stop("raw row identity/provenance failed.")
    }
    expected_methods <- c(.inference_pilot_methods(),
      paste0(.inference_pilot_methods(), "_ate"),
      "one_round_crossfit_ate_armwise",
      "one_round_crossfit_ate_hard_threshold")
    if (any(vapply(.inference_pilot_rhos(), function(rho) {
      methods <- rows$method[rows$rho == rho]
      length(methods) != 14L || anyDuplicated(methods) ||
        !setequal(methods, expected_methods)
    }, logical(1L)))) stop("raw 84-row method/rho schema failed.")
    entry <- saved$group_result$artifacts[["0"]]
    reference <- entry$direct_tate_results$one_round_crossfit
    target <- entry$data_split$t
    .ritcv_integer_scalar(reference$n_folds, "reference n_folds", 5L, 5L)
    .ritcv_integer_scalar(target$n, "target n", 1000L, 1000L)
    if (is.null(reference) || reference$n_folds != 5L ||
        !identical(reference$family, "binomial") ||
        !identical(reference$nuisance_lambda_rule, "min") ||
        !is.list(target) || target$n != 1000L ||
        any(vapply(c("A", "Y"), function(name) {
          !is.numeric(target[[name]]) || length(target[[name]]) != target$n ||
            any(!is.finite(target[[name]])) ||
            any(!target[[name]] %in% c(0, 1))
        }, logical(1L))) ||
        any(vapply(c("Z_site", "W_outcome"), function(name) {
          !is.matrix(target[[name]]) || nrow(target[[name]]) != target$n ||
            any(!is.finite(target[[name]]))
        }, logical(1L))) || !is.null(target$cv_group_id)) {
      stop("rho-zero reference/target preflight contract failed.")
    }
    expected_rho_keys <- vapply(.inference_pilot_rhos(), format, character(1L),
                                scientific = FALSE, trim = TRUE)
    if (!identical(names(saved$group_result$artifacts), expected_rho_keys)) {
      stop("saved artifact rho schema failed.")
    }
    folds <- lapply(reference$intermediates$fold_info, function(info) {
      list(original_idx = info$target_idx, n = length(info$target_idx))
    })
    indices <- lapply(folds, `[[`, "original_idx")
    if (length(folds) != 5L || any(vapply(indices, function(index) {
      !is.numeric(index) || !length(index) || any(!is.finite(index)) ||
        any(index != floor(index)) || any(index < 1L) || any(index > target$n)
    }, logical(1L))) ||
        !identical(sort(as.integer(unlist(indices))), seq_len(target$n))) {
      stop("rho-zero target outer folds do not form the required partition.")
    }
    attr(folds, ".data_ref") <- target
    for (arm in c(1L, 0L)) for (k1 in seq_len(5L)) {
      fold_result <- reference$arm_results[[paste0("mu", arm)]]$fold_results[[k1]]
      candidates <- c(list(outer = fold_result$target_only),
        fold_result$target_only_inner)
      expected_names <- c("outer", paste0("k2_", setdiff(seq_len(5L), k1)))
      if (!setequal(names(candidates), expected_names) || length(candidates) != 5L ||
          any(vapply(candidates, function(candidate) {
            !is.list(candidate) || !is.numeric(candidate$prop_scores) ||
              !is.numeric(candidate$m_pred) || any(!is.finite(candidate$prop_scores)) ||
              any(!is.finite(candidate$m_pred)) ||
              length(candidate$prop_scores) != length(candidate$m_pred)
          }, logical(1L)))) {
        stop("saved target prediction preflight schema failed.")
      }
      for (name in names(candidates)) {
        eval_fold <- if (name == "outer") k1 else as.integer(sub("^k2_", "", name))
        if (length(candidates[[name]]$prop_scores) != length(indices[[eval_fold]])) {
          stop("saved target prediction/evaluation-fold length mismatch.")
        }
      }
    }
    glm_spec <- RoCE:::resolve_glm_family(reference$family)
    diagnostics <- warning_fits <- fold_assignments <- list(); index <- 0L

    run_call <- function(model, arm, k1, k2, training, evaluation, expected,
                         shared_ps_expected = NULL) {
      state$current <- list(model = model, arm = arm, k1 = k1,
                            k2 = if (is.null(k2)) NA_integer_ else k2)
      x <- if (model == "PS") training$Z_site else
        training$W_outcome[training$A == arm, , drop = FALSE]
      y <- if (model == "PS") training$A else training$Y[training$A == arm]
      family <- if (model == "PS") "binomial" else glm_spec$glmnet_family
      .ritcv_require_standard_cv(y, model, family)
      n_cv_folds <- RoCE:::get_cv_fold_count(nrow(x), min_per_fold = 10L)
      seed <- RoCE:::.target_only_cv_seed(
        k1, k2, if (model == "PS") "propensity" else "outcome",
        if (model == "PS") 0L else arm
      )
      warnings <- character()
      cv_fit <- RoCE:::with_seed(seed, withCallingHandlers({
        prepared <- .ritcv_prepare_glmnet_args(
          x, y, family, n_cv_folds, 100L, RoCE:::GLMNET_MAX_ITER,
          NULL, paste("estimate_complement_fold_aipw", model),
          RoCE:::.make_nuisance_cv_fold_id, keep = TRUE
        )
        do.call(glmnet::cv.glmnet, prepared$args)
      }, warning = function(w) {
        message <- conditionMessage(w)
        warnings <<- c(warnings, message)
        event <- data.frame(
          event_index = length(state$warning_events) + 1L,
          model = model, arm = arm, k1 = k1,
          k2 = if (is.null(k2)) NA_integer_ else k2,
          message = message,
          condition_class = paste(class(w), collapse = ";"),
          condition_call = paste(deparse(conditionCall(w)), collapse = " "),
          stringsAsFactors = FALSE
        )
        state$warning_events[[length(state$warning_events) + 1L]] <- event
        tryInvokeRestart("muffleWarning")
      }))
      prediction <- as.numeric(predict(
        cv_fit,
        newx = if (model == "PS") evaluation$Z_site else evaluation$W_outcome,
        s = "lambda.min", type = "response"
      ))
      prediction <- if (model == "PS") RoCE:::clip_propensity(prediction) else
        RoCE:::clip_outcome_pred(prediction, reference$family)
      prediction_error <- .ritcv_validate_prediction(prediction, expected)
      shared_ps_error <- if (is.null(shared_ps_expected)) NA_real_ else
        .ritcv_validate_prediction(prediction, shared_ps_expected)
      actual_fold_id <- as.integer(cv_fit$foldid)
      if (length(actual_fold_id) != nrow(x) || anyNA(actual_fold_id) ||
          any(!actual_fold_id %in% seq_len(n_cv_folds))) {
        stop("cv.glmnet did not retain a valid realized fold assignment.")
      }
      result <- data.frame(
        model = model, arm = arm, k1 = k1,
        k2 = if (is.null(k2)) NA_integer_ else k2,
        seed = seed, n_train = nrow(x), n_cv_folds = n_cv_folds,
        fold_assignment_mode = "internal_random_nfolds",
        lambda_path_length = length(cv_fit$lambda),
        lambda_path_min = min(cv_fit$lambda), lambda_path_max = max(cv_fit$lambda),
        lambda_min = cv_fit$lambda.min, lambda_1se = cv_fit$lambda.1se,
        selected_lambda = cv_fit$lambda.min,
        full_model_jerr = if (is.null(cv_fit$glmnet.fit$jerr)) NA_integer_ else
          as.integer(cv_fit$glmnet.fit$jerr),
        fit_preval_rows = nrow(cv_fit$fit.preval),
        fit_preval_columns = ncol(cv_fit$fit.preval),
        prediction_length = length(prediction), saved_prediction_length = length(expected),
        prediction_max_error = prediction_error,
        shared_ps_prediction_max_error = shared_ps_error,
        warning_count = length(warnings),
        warning = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE
      )
      list(
        diagnostic = result,
        fold_id = actual_fold_id,
        warning_fit = if (length(warnings)) cv_fit else NULL
      )
    }

    for (arm in c(1L, 0L)) for (k1 in seq_len(reference$n_folds)) {
      for (k2_value in c(NA_integer_, setdiff(seq_len(reference$n_folds), k1))) {
        k2 <- if (is.na(k2_value)) NULL else k2_value
        excluded <- if (is.null(k2)) k1 else c(k1, k2)
        training <- RoCE:::combine_folds(folds, setdiff(seq_len(reference$n_folds), excluded))
        evaluation <- RoCE:::materialize_fold(folds, if (is.null(k2)) k1 else k2)
        stored <- if (is.null(k2)) reference$arm_results[[paste0("mu", arm)]]$
          fold_results[[k1]]$target_only else reference$arm_results[[paste0("mu", arm)]]$
          fold_results[[k1]]$target_only_inner[[paste0("k2_", k2)]]
        models <- if (arm == 1L) c("PS", "OR") else "OR"
        for (model in models) {
          index <- index + 1L
          replayed <- run_call(model, arm, k1, k2, training, evaluation,
            if (model == "PS") stored$prop_scores else stored$m_pred,
            if (model == "PS") {
              other <- if (is.null(k2)) reference$arm_results$mu0$
                fold_results[[k1]]$target_only else reference$arm_results$mu0$
                fold_results[[k1]]$target_only_inner[[paste0("k2_", k2)]]
              other$prop_scores
            } else NULL)
          diagnostics[[index]] <- replayed$diagnostic
          state$diagnostics[[length(state$diagnostics) + 1L]] <-
            replayed$diagnostic
          call_key <- paste(
            model, arm, k1, if (is.null(k2)) "NA" else k2, sep = "_"
          )
          fold_assignments[[call_key]] <- replayed$fold_id
          if (!is.null(replayed$warning_fit)) {
            warning_fits[[call_key]] <- replayed$warning_fit
          }
        }
      }
    }
    diagnostics <- do.call(rbind, diagnostics)
    if (nrow(diagnostics) != 75L || sum(diagnostics$model == "PS") != 25L ||
        sum(diagnostics$model == "OR") != 50L ||
        any(diagnostics$prediction_max_error > 1e-10) ||
        any(diagnostics$shared_ps_prediction_max_error[
          diagnostics$model == "PS"
        ] > 1e-10) ||
        any(!is.na(diagnostics$shared_ps_prediction_max_error[
          diagnostics$model == "OR"
        ]))) {
      stop("75-call replay completeness/identity gate failed.")
    }
    warning_events <- if (length(state$warning_events))
      do.call(rbind, state$warning_events) else .ritcv_empty_warning_events()
    list(diagnostics = diagnostics, warning_events = warning_events,
         warning_fits = warning_fits,
         fold_assignments = fold_assignments,
         task_id = task_id, sim_id = sim_id, bundle_hash = bundle_hash,
         audit = audit, audit_hash = audit_hash,
         package_hash = actual_package)
  }, error = identity)

  if (inherits(attempt, "error")) {
    .ritcv_publish_failure(output, attempt, state, roce_sha256_file,
                           roce_write_atomic_directory)
    stop("independent target CV replay failed; failure artifact published: ",
         conditionMessage(attempt), call. = FALSE)
  }
  roce_write_atomic_directory(output, function(stage) {
    if (.ritcv_check_sha(bundle, c(
      "results.csv", "inference_audit.csv", "diagnostic_qc.csv",
      "artifacts.rds", "metadata.txt"
    ), roce_sha256_file) != attempt$bundle_hash ||
        .ritcv_check_sha(attempt$audit, c(
          "audit_passed.txt", "summary.csv", "numerical_checks.csv",
          "nuisance_diagnostics.csv"
        ), roce_sha256_file) != attempt$audit_hash ||
        .weight_calibration_installed_package_fingerprint(
          package_library, roce_sha256_file
        ) != attempt$package_hash ||
        .weight_calibration_files_fingerprint(
          import_files, roce_sha256_file
        ) != import_hash ||
        roce_sha256_file(
          "diagnosis/tate_common_weight/replay_independent_target_cv.R"
        ) != script_hash) stop("replay inputs changed before publication.")
    write.csv(attempt$diagnostics, file.path(stage, "target_cv_replay.csv"), row.names = FALSE)
    write.csv(attempt$warning_events, file.path(stage, "warning_events.csv"),
              row.names = FALSE)
    saveRDS(attempt$fold_assignments, file.path(stage, "cv_fold_assignments.rds"))
    saveRDS(attempt$warning_fits, file.path(stage, "warning_cv_fits.rds"))
    writeLines(c("independent_target_cv_replay=passed", "replay_calls=75",
      "saved_partition_identity_verified=FALSE",
      "warning_lambda_index_is_master_grid_index=FALSE",
      paste0("task_id=", attempt$task_id), paste0("sim_id=", attempt$sim_id),
      paste0("bundle_checksum_manifest_fingerprint=", attempt$bundle_hash),
      paste0("seed_audit_checksum_manifest_fingerprint=", attempt$audit_hash),
      paste0("package_fingerprint=", attempt$package_hash),
      paste0("replay_script_fingerprint=", script_hash),
      paste0("import_dependency_fingerprint=", import_hash),
      paste0("R_version=", R.version.string),
      paste0("RoCE_version=", as.character(utils::packageVersion("RoCE"))),
      paste0("glmnet_version=", as.character(utils::packageVersion("glmnet")))),
      file.path(stage, "metadata.txt"))
    files <- c("target_cv_replay.csv", "warning_events.csv",
               "cv_fold_assignments.rds",
               "warning_cv_fits.rds", "metadata.txt")
    writeLines(paste(vapply(file.path(stage, files), roce_sha256_file, character(1L)),
                     files, sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "independent target CV replay")
  invisible(attempt$diagnostics)
}

if (sys.nframe() == 0L) replay_independent_target_cv_main()

#!/usr/bin/env Rscript

# Standalone, read-only implementation audit for one independently generated
# inference-pilot seed. This file is deliberately outside the frozen runner's
# workflow fingerprint.

.iia_metadata <- function(path) {
  x <- readLines(path, warn = FALSE); p <- regexpr("=", x, fixed = TRUE)
  if (!length(x) || any(p < 2L)) stop("malformed metadata", call. = FALSE)
  k <- substr(x, 1L, p - 1L); v <- substr(x, p + 1L, nchar(x))
  if (anyDuplicated(k)) stop("duplicate metadata keys", call. = FALSE)
  stats::setNames(v, k)
}

.iia_manifest_path <- function(bundle) {
  d <- normalizePath(bundle, mustWork = TRUE)
  repeat {
    candidate <- file.path(d, "manifest.csv")
    if (file.exists(candidate)) return(normalizePath(candidate))
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("no ancestor manifest.csv found", call. = FALSE)
}

.iia_strict_integer <- function(x, name, lower = 1L, upper = .Machine$integer.max) {
  if (!is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      is.logical(x) || is.complex(x) ||
      !typeof(x) %in% c("integer", "double", "character")) {
    stop(name, " must be one strict integer", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(x))
  if (length(value)!=1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper)
    stop(name," must be one strict integer",call.=FALSE)
  as.integer(value)
}

.iia_validate_saved_task <- function(observed, expected) {
  if (!is.data.frame(observed) || nrow(observed) != 1L ||
      !identical(names(observed), names(expected)) ||
      any(vapply(names(expected), function(name) {
        value <- observed[[name]]
        !is.atomic(value) || !is.null(dim(value)) || is.object(value) ||
          anyNA(value) ||
          !identical(as.character(value), as.character(expected[[name]]))
      }, logical(1L)))) {
    stop("saved task disagrees with canonical manifest task", call. = FALSE)
  }
  invisible(TRUE)
}

.iia_binary_vector <- function(x, n) {
  (is.numeric(x) || is.logical(x)) && is.null(dim(x)) && !is.object(x) &&
    length(x) == n && all(is.finite(x)) && all(x %in% c(0, 1))
}

.iia_numeric_errors <- function(values, expected_names) {
  if (!is.numeric(values) || length(values)!=length(expected_names) ||
      !identical(names(values),expected_names) || any(!is.finite(values)) ||
      any(values < 0)) stop("numerical checks must be fixed, finite, and nonnegative",call.=FALSE)
  values
}

.iia_validate_data_fit <- function(entry, fit) {
  sites <- c("t","s1","s2"); data <- entry$data_split
  if (!is.list(data) || !identical(names(data),sites)) stop("invalid data_split sites",call.=FALSE)
  for (site in sites) {
    x <- data[[site]]
    if (!is.list(x) || !all(c("n","X","A","Y")%in%names(x)) ||
        .iia_strict_integer(x$n,paste(site,"n"))!=1000L ||
        !is.matrix(x$X) || !is.numeric(x$X) ||
        !identical(dim(x$X),c(1000L,100L)) || any(!is.finite(x$X)) ||
        !.iia_binary_vector(x$A, 1000L) || !.iia_binary_vector(x$Y, 1000L))
      stop("invalid data_split site schema: ",site,call.=FALSE)
  }
  if (!is.numeric(fit$variance) || length(fit$variance)!=1L ||
      !is.null(dim(fit$variance)) || !is.null(dim(fit$all_phi_agg)) ||
      !is.finite(fit$variance) || fit$variance<=0 || !is.numeric(fit$all_phi_agg) ||
      length(fit$all_phi_agg)!=3000L || any(!is.finite(fit$all_phi_agg)) ||
      .iia_strict_integer(fit$N_all, "fit N_all") != 3000L ||
      .iia_strict_integer(fit$n_folds, "fit n_folds") != 5L ||
      .iia_strict_integer(fit$n_sites, "fit n_sites") != 3L)
    stop("fit variance/influence/sample-size schema invalid",call.=FALSE)
  invisible(TRUE)
}

.iia_verify_manifest <- function(directory, expected, sha256_file) {
  entries <- list.files(directory, all.files = TRUE, no.. = TRUE,
                        recursive = FALSE, include.dirs = TRUE)
  if (!setequal(entries, c(expected, "sha256.txt")) ||
      length(entries) != length(expected) + 1L)
    stop("bundle directory has extra or missing payload entries", call. = FALSE)
  lines <- readLines(file.path(directory, "sha256.txt"), warn = FALSE)
  hashes <- substr(lines, 1L, 64L); files <- substring(lines, 67L)
  if (length(lines) != length(expected) || anyDuplicated(files) ||
      !setequal(files, expected) || any(substr(lines, 65L, 66L) != "  ") ||
      any(!grepl("^[0-9a-f]{64}$", hashes)))
    stop("bundle checksum manifest has the wrong exact payload", call. = FALSE)
  observed <- vapply(file.path(directory, files), sha256_file, character(1L))
  if (!identical(unname(observed), hashes)) stop("bundle checksum mismatch", call. = FALSE)
  invisible(stats::setNames(hashes, files))
}

.iia_classify_bundle <- function(bundle, sha256_file) {
  success <- c("results.csv", "inference_audit.csv", "diagnostic_qc.csv",
               "artifacts.rds", "metadata.txt")
  if (file.exists(file.path(bundle, "attempt_failure.csv"))) {
    .iia_verify_manifest(bundle, c("attempt_failure.csv", "metadata.txt"),
                         sha256_file)
    failure <- read.csv(file.path(bundle, "attempt_failure.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("status", "failure_stage", "failure_message")
    if (nrow(failure) != 1L || !all(required %in% names(failure)) ||
        failure$status[[1L]] != "failed" ||
        is.na(failure$failure_stage[[1L]]) ||
        !nzchar(failure$failure_stage[[1L]]) ||
        is.na(failure$failure_message[[1L]]) ||
        !nzchar(failure$failure_message[[1L]]))
      stop("failed-attempt record lacks one explicit stage/reason", call. = FALSE)
    return(list(status = "failed", stage = failure$failure_stage[[1L]],
                reason = failure$failure_message[[1L]]))
  }
  .iia_verify_manifest(bundle, success, sha256_file)
  list(status = "complete", stage = NULL, reason = NULL)
}

.iia_nuisance <- function(fit, rho) {
  d <- fit$nuisance_fit_diagnostics
  if (!is.data.frame(d) || !nrow(d)) stop("missing nuisance diagnostics", call. = FALSE)
  fields <- c("initial_outcome_degenerate", "target_only_outcome_degenerate",
    "initial_dr_nonconverged", "calibrated_dr_nonconverged",
    "calibrated_outcome_nonconverged", "initial_dr_line_search_failures",
    "calibrated_dr_line_search_failures", "calibrated_outcome_line_search_failures",
    "initial_dr_support_floor_applied", "calibrated_dr_support_floor_applied")
  missing <- setdiff(fields, names(d))
  if (length(missing)) stop("nuisance diagnostics missing hard fields", call. = FALSE)
  if (any(vapply(d[fields],function(x)!is.numeric(x) || any(!is.finite(x)) ||
      any(x<0) || any(x!=floor(x)),logical(1L))))
    stop("hard nuisance fields must be known finite nonnegative integers",call.=FALSE)
  data.frame(rho = rho, rows = nrow(d), field = fields,
    total = vapply(d[fields], sum, numeric(1L)),
    stringsAsFactors = FALSE)
}

audit_independent_inference_pilot_main <- function(args = commandArgs(trailingOnly=TRUE)) {
  if (length(args) != 3L) stop(paste("usage: audit_independent_inference_pilot.R",
    "BUNDLE_DIR PACKAGE_LIBRARY OUTPUT_DIR"), call. = FALSE)
  bundle <- normalizePath(args[[1L]], mustWork = TRUE)
  package_library <- normalizePath(args[[2L]], mustWork = TRUE)
  output <- file.path(normalizePath(dirname(args[[3L]]), mustWork = FALSE), basename(args[[3L]]))
  if (file.exists(output) || dir.exists(output)) stop("audit output already exists", call. = FALSE)
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
  source("diagnosis/tate_common_weight/audit_full_refit_aggregation_intermediates.R")
  .libPaths(c(package_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  if (!identical(dirname(normalizePath(find.package("RoCE"))), package_library))
    stop("RoCE not loaded from PACKAGE_LIBRARY", call. = FALSE)

  success_files <- c("results.csv","inference_audit.csv","diagnostic_qc.csv",
                     "artifacts.rds","metadata.txt")
  classification <- .iia_classify_bundle(bundle, roce_sha256_file)
  if (classification$status == "failed")
    stop("refusing success audit for failed attempt: stage=", classification$stage,
         "; reason=", classification$reason, call. = FALSE)
  bundle_manifest_fp <- roce_sha256_file(file.path(bundle,"sha256.txt"))
  metadata <- .iia_metadata(file.path(bundle, "metadata.txt"))
  manifest_path <- .iia_manifest_path(bundle)
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  task_id <- .iia_strict_integer(metadata[["task_id"]],"metadata task_id",1L,100L)
  task <- roce_validate_inference_pilot_manifest(manifest, task_id)
  expected <- c(package_fingerprint="2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    workflow_fingerprint="9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4")
  if (metadata[["independent_inference_pilot"]] != "complete" ||
      metadata[["package_fingerprint"]] != expected[[1L]] ||
      metadata[["workflow_fingerprint"]] != expected[[2L]] ||
      metadata[["manifest_fingerprint"]] != roce_sha256_file(manifest_path) ||
      .iia_strict_integer(metadata[["sim_id"]],"metadata sim_id",10001L,10100L) != task$sim_id ||
      metadata[["n_weight_bootstrap"]] != "1000")
    stop("pilot metadata/provenance gate failed", call. = FALSE)
  installed <- .weight_calibration_installed_package_fingerprint(package_library, roce_sha256_file)
  if (installed != expected[[1L]]) stop("installed package fingerprint mismatch", call. = FALSE)
  package_gate <- readLines(file.path(package_library,"audit_tests_passed.txt"),warn=FALSE)
  if (!"package_tests=passed" %in% package_gate ||
      !paste0("package_fingerprint=",installed) %in% package_gate ||
      !"package_source_fingerprint=8e17115d87a5cfa0ba9f4051c4018989eb9c3b96067d5f9eeece792edd4b82ac" %in% package_gate)
    stop("installed package source/test gate is absent or inconsistent",call.=FALSE)
  workflow_files <- c("diagnosis/tate_common_weight/run_independent_inference_pilot.R",
    "diagnosis/tate_common_weight/run_independent_inference_pilot.sh",
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R",
    "scripts/slurm/atomic_output.R","scripts/slurm/result_provenance.R",
    "scripts/slurm/resource_topology.R","scripts/slurm/direct_tate_task_helpers.R",
    "scripts/slurm/simulation_qc_policy.R","scripts/slurm/package_library_utils.sh")
  if (.weight_calibration_files_fingerprint(workflow_files,roce_sha256_file)!=expected[[2L]])
    stop("live nine-file workflow fingerprint mismatch",call.=FALSE)

  results <- read.csv(file.path(bundle,"results.csv"),check.names=FALSE,na.strings=c("","NA"))
  inference <- read.csv(file.path(bundle,"inference_audit.csv"),check.names=FALSE)
  qc <- read.csv(file.path(bundle,"diagnostic_qc.csv"),check.names=FALSE)
  saved <- readRDS(file.path(bundle,"artifacts.rds"))
  .iia_validate_saved_task(saved$task, task)
  locked_args <- list(sim_id=as.integer(task$sim_id),n_total=3000L,K=2L,p=100L,
    config="C1",methods=.inference_pilot_methods(),nlambda_init=100L,nuisance_lambda_rule="min",
    estimand_type="superpopulation",outcome_type="binary",n_folds=5L,
    aggregation_lambda=1,n_bootstrap=5000L,n_weight_bootstrap=1000L,
    M_tau=5,M_tau_inference=5,estimate_ate=TRUE,
    include_hard_threshold_diagnostic=TRUE,dgp_type="face")
  if (!is.list(saved$simulation_args) || any(vapply(names(locked_args),function(n)
      !identical(saved$simulation_args[[n]],locked_args[[n]]),logical(1L))))
    stop("saved simulation_args disagree with frozen scientific settings",call.=FALSE)
  expected_methods <- c(.inference_pilot_methods(),paste0(.inference_pilot_methods(),"_ate"),
    "one_round_crossfit_ate_armwise","one_round_crossfit_ate_hard_threshold")
  if (nrow(results)!=84L || nrow(inference)!=6L || nrow(qc)!=84L ||
      any(results$sim_id != task$sim_id) || any(results$task_id != task$task_id) ||
      !identical(sort(unique(results$rho)),.inference_pilot_rhos()) ||
      any(vapply(.inference_pilot_rhos(),function(rho)
        !setequal(results$method[results$rho==rho],expected_methods),logical(1L))) ||
      any(results$package_fingerprint != expected[[1L]]) ||
      any(results$workflow_fingerprint != expected[[2L]]) ||
      any(results$manifest_fingerprint != metadata[["manifest_fingerprint"]]) ||
      any(qc$implementation_failed))
    stop("result/QC shape or implementation gate failed", call. = FALSE)
  rebuilt <- roce_build_inference_audit(results)
  if (!isTRUE(all.equal(inference, rebuilt, tolerance=1e-12, check.attributes=FALSE)))
    stop("derived inference audit does not rebuild", call. = FALSE)
  keys <- vapply(.inference_pilot_rhos(), format, character(1L), scientific=FALSE,trim=TRUE)
  if (!identical(names(saved$group_result$artifacts),keys)) stop("artifact rho keys differ",call.=FALSE)

  checks <- list(); nuis <- list(); index <- 0L
  for (rho in .inference_pilot_rhos()) {
    key <- format(rho,scientific=FALSE,trim=TRUE)
    entry <- saved$group_result$artifacts[[key]]
    fit <- entry$direct_tate_results$one_round_crossfit
    RoCE:::.validate_weight_bootstrap_result(fit)
    .iia_validate_data_fit(entry,fit)
    row <- results[results$rho==rho & results$method=="one_round_crossfit_ate",]
    fold <- do.call(rbind,lapply(seq_len(fit$n_folds),function(k).afri_manual_fold(fit,key,k)))
    arm <- .afri_arm_audit(fit,key)
    sizes <- c(fit$intermediates$sample_sizes$n_t,fit$intermediates$sample_sizes$n_source)
    values <- c(result_estimate=abs(row$estimate-fit$estimate),result_se=abs(row$se-fit$se),
      mean_phi=abs(fit$estimate-mean(fit$all_phi_agg)),
      site_variance=abs(fit$variance-RoCE:::.multisite_pseudovalue_variance(fit$all_phi_agg,sizes)),
      arm_contrast=max(c(arm$fold_arm_contrast_error,arm$inner_arm_contrast_error)),
      fold_aggregation=max(arm$fold_aggregation_identity_error),
      wald=max(fold$wald_reconstruction_error),penalty=max(fold$penalty_reconstruction_error),
      variance_api=max(fold$variance_api_reconstruction_error),
      psd_ridge=max(fold$psd_ridge_reconstruction_error),
      optimizer=max(fold$optimizer_rerun_weight_error),kkt=max(fold$normalized_kkt_residual_max))
    values <- .iia_numeric_errors(values,c("result_estimate","result_se","mean_phi",
      "site_variance","arm_contrast","fold_aggregation","wald","penalty",
      "variance_api","psd_ridge","optimizer","kkt"))
    index <- index+1L
    checks[[index]] <- data.frame(rho=rho,check=names(values),maximum_error=values,row.names=NULL)
    nuis[[index]] <- .iia_nuisance(fit,rho)
  }
  numerical <- do.call(rbind,checks); nuisance <- do.call(rbind,nuis)
  if (any(numerical$maximum_error[numerical$check!="kkt"]>1e-10) ||
      any(numerical$maximum_error[numerical$check=="kkt"]>1e-5) ||
      any(nuisance$total[nuisance$field %in% c("initial_dr_nonconverged",
        "calibrated_dr_nonconverged","calibrated_outcome_nonconverged",
        "initial_dr_line_search_failures","calibrated_dr_line_search_failures",
        "calibrated_outcome_line_search_failures")] > 0))
    stop("numerical or hard nuisance audit failed", call. = FALSE)
  audit_driver_fp <- roce_sha256_file("diagnosis/tate_common_weight/audit_independent_inference_pilot.R")
  afri_fp <- roce_sha256_file("diagnosis/tate_common_weight/audit_full_refit_aggregation_intermediates.R")
  runner_fp <- roce_sha256_file("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
  calibration_fp <- roce_sha256_file("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  atomic_fp <- roce_sha256_file("scripts/slurm/atomic_output.R")
  provenance_fp <- roce_sha256_file("scripts/slurm/result_provenance.R")
  slurm_ids <- unique(as.character(results$slurm_job_id))
  summary <- data.frame(task_id=task$task_id,sim_id=task$sim_id,n_rho=6L,
    result_rows=nrow(results),max_numerical_error=max(numerical$maximum_error),
    implementation_failures=sum(qc$implementation_failed),
    coverage_status="diagnostic_only_not_a_coverage_gate",
    warning_capture_complete=FALSE,warning_count=NA_integer_,
    runtime_slurm_job_ids=paste(slurm_ids,collapse=";"),stringsAsFactors=FALSE)
  roce_write_atomic_directory(output,function(stage){
    .iia_verify_manifest(bundle,success_files,roce_sha256_file)
    if (roce_sha256_file(file.path(bundle,"sha256.txt")) != bundle_manifest_fp)
      stop("bundle checksum manifest changed before audit publication",call.=FALSE)
    write.csv(summary,file.path(stage,"summary.csv"),row.names=FALSE)
    write.csv(numerical,file.path(stage,"numerical_checks.csv"),row.names=FALSE)
    write.csv(nuisance,file.path(stage,"nuisance_diagnostics.csv"),row.names=FALSE)
    writeLines(c("independent_inference_seed_audit=passed",paste0("bundle_directory=",bundle),
      paste0("package_fingerprint=",metadata[["package_fingerprint"]]),
      paste0("workflow_fingerprint=",metadata[["workflow_fingerprint"]]),
      paste0("manifest_fingerprint=",metadata[["manifest_fingerprint"]]),
      paste0("bundle_checksum_manifest_fingerprint=",bundle_manifest_fp),
      paste0("audit_driver_fingerprint=",audit_driver_fp),
      paste0("afri_dependency_fingerprint=",afri_fp),
      paste0("runner_dependency_fingerprint=",runner_fp),
      paste0("calibration_dependency_fingerprint=",calibration_fp),
      paste0("atomic_dependency_fingerprint=",atomic_fp),
      paste0("provenance_dependency_fingerprint=",provenance_fp),
      paste0("task_id=",task$task_id),paste0("sim_id=",task$sim_id),
      "coverage_status=diagnostic_only_not_passed",
      "warning_capture_complete=FALSE","warning_count=NA",
      paste0("runtime_slurm_job_ids=",paste(slurm_ids,collapse=";")),
      "warning_capture_status=runner_has_no_complete_warning field; absence cannot establish zero warnings"),
      file.path(stage,"audit_passed.txt"))
    files<-c("audit_passed.txt","summary.csv","numerical_checks.csv","nuisance_diagnostics.csv")
    h<-vapply(file.path(stage,files),roce_sha256_file,character(1L))
    writeLines(paste(h,files,sep="  "),file.path(stage,"sha256.txt"))
  },caller="independent inference seed audit")
  message("[passed] independent inference seed audit: ",output)
}

if (sys.nframe()==0L) audit_independent_inference_pilot_main()

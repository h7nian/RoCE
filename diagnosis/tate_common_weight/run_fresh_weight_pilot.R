#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: run_fresh_weight_pilot.R LIBRARY SIM_ID OUTPUT")
  sim_id<-suppressWarnings(as.integer(args[2]));stopifnot(!is.na(sim_id),sim_id %in% 20001:20005)
  if(file.exists(args[3]))stop("output already exists")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("scripts/slurm/direct_tate_task_helpers.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  for(file in c("quadratic_bias_weights.R","quadratic_weight_influence.R","balanced_fold_gradient.R"))
    source(file.path("diagnosis/tate_common_weight",file))
  fingerprint<-.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)
  stopifnot(fingerprint=="6a715f3a41f2f42ef58bc08310d4ecca061524b0ab555d52a6221101a63a0a93")
  release<-roce_claim_task_lock(paste0(args[3],".lock"),c(paste0("sim_id=",sim_id),paste0("job_id=",Sys.getenv("SLURM_JOB_ID"))))
  on.exit(release(),add=TRUE)
  # This API validates an unused comparison-bootstrap count even with no
  # comparison methods. Abort if its resampling implementation is ever reached.
  trace(".multiplier_bootstrap_se",where=asNamespace("RoCE"),
    tracer=quote(stop("bootstrap is disabled in the fresh analytic pilot")),print=FALSE)
  on.exit(untrace(".multiplier_bootstrap_se",where=asNamespace("RoCE")),add=TRUE)
  rhos<-c(0,.5,1,1.5,2,2.5)
  simulation_args<-list(sim_id=sim_id,n_total=3000L,K=2L,p=100L,config="C1",
    methods="one_round_crossfit",verbose=FALSE,n_cores_internal=1L,nlambda_init=100L,
    nuisance_lambda_rule="min",estimand_type="superpopulation",outcome_type="binary",
    n_folds=5L,aggregation_lambda=1,n_bootstrap=2L,n_weight_bootstrap=0L,
    M_tau=5,M_tau_inference=5,estimate_ate=TRUE,parallel_treatment_arms=TRUE,dgp_type="face")
  grouped<-NULL;gradients<-rows<-list();warnings<-character()
  failure<-tryCatch(withCallingHandlers({
    grouped<-RoCE:::.run_face_rho_group(simulation_args,rho_values=rhos,changed_sources="s1",
      artifact_rhos=rhos,positive_rho_workers=1L)
    stopifnot(identical(names(grouped$artifacts),as.character(rhos)))
    for(rho in rhos) {
      artifact<-grouped$artifacts[[as.character(rho)]];fit<-artifact$direct_tate_results$one_round_crossfit
      q<-.quadratic_weight_influence(fit);data<-artifact$data_split
      balanced<-lapply(seq_along(data),function(site) {
        folds<-integer(data[[site]]$n)
        for(k in 1:5) {
          f<-fit$intermediates$fold_info[[k]];ids<-if(site==1L)f$target_idx else f$source_idx[[site-1L]]
          folds[ids]<-k
        }
        .balanced_arm_fold_gradient(q$gradient[[site]],data[[site]]$A,folds)
      })
      diagnostics<-fit$nuisance_fit_diagnostics
      hard<-grep("nonconverged|degenerate",names(diagnostics),value=TRUE)
      stopifnot(all(as.matrix(diagnostics[hard])==0),abs(mean(fit$all_phi_tau)-fit$estimate)<1e-10)
      estimates<-c(target_only=fit$target_only$estimate,original_common_weights=fit$estimate,
        quadratic_fixed_weights=q$estimate,quadratic_analytic_weights=q$estimate)
      standard_errors<-c(fit$target_only$se,fit$se,sqrt(q$fixed_variance),
        sqrt(sum(vapply(balanced,`[[`,numeric(1L),"variance"))))
      truth<-artifact$tate_truth;lower<-estimates-qnorm(.975)*standard_errors;upper<-estimates+qnorm(.975)*standard_errors
      rows[[as.character(rho)]]<-data.frame(sim_id,rho,method=names(estimates),estimate=unname(estimates),
        se=standard_errors,truth,bias=unname(estimates)-truth,ci_lower=unname(lower),ci_upper=unname(upper),
        coverage=unname(lower<=truth & upper>=truth),inference_validated=FALSE)
      gradients[[as.character(rho)]]<-list(weight_layer=q,balanced=balanced)
    }
    NULL
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
  roce_write_atomic_directory(args[3],function(stage) {
    saveRDS(list(grouped=grouped,gradients=gradients,failure=failure,warnings=warnings),file.path(stage,"artifacts.rds"))
    if(is.null(failure)) {
      results<-do.call(rbind,rows);stopifnot(nrow(results)==24L)
      write.csv(results,file.path(stage,"results.csv"),row.names=FALSE)
    }
    writeLines(c(paste0("sim_id=",sim_id),paste0("package_fingerprint=",fingerprint),
      paste0("failed=",!is.null(failure)),"experiment=fresh_analytic_weight_pilot",
      "bootstrap_draws=0","comparison_methods_requested=FALSE","nuisance_refits_per_replication=TRUE",
      "quadratic_power=0.75","warning_capture_complete=FALSE","inference_validated=FALSE",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/run_fresh_weight_pilot.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="fresh independent analytic weight pilot")
  if(!is.null(failure))stop(conditionMessage(failure))
}
if(sys.nframe()==0L)main()

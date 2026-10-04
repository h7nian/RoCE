#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: integration_root new_scratch_output")
root <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new scratch output")
library_path <- file.path(root,"Rlib")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
dir.create(output,recursive=TRUE)
configuration <- list(n_sims=2L,n_total_vec=2000L,K_vec=1L,p_vec=4L,configs="C2",
  n_target=1000L,n_source_sizes=1000L,n_folds=4L,nlambda_init=100L,
  dgp_type="bounded",methods="one_round_crossfit",estimate_ate=TRUE,
  target_nuisance_method="hou_calibrated",source_validation_method="calibrated",
  crossfit_layers=2L,calibration_layout="compact",nuisance_solver="proximal_newton",
  nuisance_tol=1e-10,M_tau=12,M_tau_inference=12,
  calibration_control=list(recipe="score_derivative",target_propensity_initialization="calibrated",target_radius=12),
  additional_aggregation_modes=c("separate_arms","joint_tate"),checkpoint_config=NULL,
  nested_parallel=FALSE,parallel_strategy="outer_only",verbose_every=100L)
previous <- options(RoCE.nuisance_cv_certificate=NULL)
reference <- do.call(run_simulation_study,c(configuration,list(n_cores=1L,nuisance_cv_certificate=FALSE)))
serial <- do.call(run_simulation_study,c(configuration,list(n_cores=1L,nuisance_cv_certificate=TRUE)))
socket <- do.call(run_simulation_study,c(configuration,list(n_cores=2L,nuisance_cv_certificate=TRUE)))
for(result in list(reference,serial,socket)) {
  stopifnot(nrow(result)>0L,"nuisance_cv_certificate"%in%names(result),
    length(result$nuisance_cv_certificate)==nrow(result),!anyNA(result$nuisance_cv_certificate))
}
stopifnot(is.null(getOption("RoCE.nuisance_cv_certificate")),!any(reference$nuisance_cv_certificate),
  all(serial$nuisance_cv_certificate),all(socket$nuisance_cv_certificate))
standardize <- function(rows) {
  rows <- rows[order(rows$sim_id,rows$method,rows$estimand_scope),
    c("sim_id","method","estimand_scope","estimate","se","truth")]
  rownames(rows) <- NULL
  rows
}
stopifnot(isTRUE(all.equal(standardize(reference),standardize(serial),tolerance=1e-10)),
  isTRUE(all.equal(standardize(reference),standardize(socket),tolerance=1e-10)))
saveRDS(list(configuration=configuration,reference=reference,serial=serial,socket=socket),
  file.path(output,"checked_results.rds"))
writeLines(c("SOCKET_CERTIFICATE_PROPAGATION_AND_PARITY_PASSED",
  "C2,p4,K1,1000/site,4 folds,two layers,100 lambdas; two repeats,three aggregation modes.",
  "Reference off versus serial on and two-worker PSOCK on; effective metadata and restored parent scope checked."),
  file.path(output,"checks.txt"))
options(previous)

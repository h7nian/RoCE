#!/usr/bin/env Rscript
# Complete path and refit comparison; no production campaign uses this flag.
args <- commandArgs(trailingOnly=TRUE)
if(!length(args)%in%c(4L,5L)) stop("Usage: checked_root reference_library workflow_directory output [reference|resume]")
check <- normalizePath(args[1L],mustWork=TRUE)
reference_library <- normalizePath(args[2L],mustWork=TRUE)
workflow <- normalizePath(args[3L],mustWork=TRUE)
output <- args[4L]
reference_mode <- length(args)==5L && args[5L]=="reference"
resume_mode <- length(args)==5L && args[5L]=="resume"
if(length(args)==5L && !args[5L]%in%c("reference","resume")) stop("Unknown validation mode")
library_path <- if(reference_mode) reference_library else file.path(check,"Rlib")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")

refit <- function(specification,lambda) {
  a <- specification$arguments
  if(specification$type=="initial") {
    RoCE::fit_initial_density_ratio(a$Z_site,a$A_source,a$mean_phi,lambda=lambda,
      A_val=a$A_val,M_tau=12,tol=1e-10)
  } else {
    RoCE::fit_unified_density_ratio(a$Z_site,a$A_source,a$mean_grad_psi,a$alpha_init,
      lambda=lambda,A_val=a$A_val,M_tau=12,tol=1e-10,
      W_outcome=specification$outcome_design,calibrated=specification$type=="calibrated",
      family_int=1L,link_int=1L,truncate_initial_outcome=FALSE)
  }
}

if(reference_mode) {
  specifications <- readRDS(file.path(output,"specifications.rds"))
  results <- lapply(specifications,function(specification) {
    Sys.setenv(ROCE_NUISANCE_SOLVER=specification$solver)
    RoCE:::set_nuisance_solver_cpp(specification$solver)
    cv <- do.call(getFromNamespace(specification$function_name,"RoCE"),specification$arguments)
    list(cv=cv,fit=refit(specification,cv$lambda_min))
  })
  saveRDS(results,file.path(output,"reference_results.rds"))
  quit(status=0L)
}

if(dir.exists(output) && !resume_mode) stop("Use a new output directory")
if(!dir.exists(output) && resume_mode) stop("No validation checkpoint to resume")
if(resume_mode) {
  specifications <- readRDS(file.path(output,"specifications.rds"))
  stopifnot(file.exists(file.path(output,"reference_results.rds")))
} else {
dir.create(output,recursive=TRUE)
specifications <- list()
cases <- data.frame(p=c(10L,200L,200L,200L),config=c("C1","C1","C2","C3"),
                    arm=c(1L,0L,1L,0L),seed=2841L+0:3)
for(case_index in seq_len(nrow(cases))) {
  case <- cases[case_index,]
  set.seed(case$seed)
  data <- split_data_by_site(generate_bounded_data(n_target=1000L,n_source_sizes=1000L,
    K=1L,p=case$p,config=case$config))
  folds <- RoCE:::build_crossfit_folds(data,10L)
  source <- RoCE:::combine_folds(folds$source_folds$s1,4:10)
  target <- RoCE:::combine_folds(folds$target_folds,4:10)
  beta <- RoCE::fit_initial_outcome(target$W_outcome,target$Y,target$A,
    A_val=case$arm,nlambda=100L)
  for(type in c("initial","refined","calibrated")) {
    outcome_design <- source$W_outcome
    target_design <- target$W_outcome
    alpha <- as.numeric(beta)
    if(type=="calibrated") {
      outcome_design <- matrix(drop(cbind(1,outcome_design)%*%alpha),ncol=1L)
      target_design <- matrix(drop(cbind(1,target_design)%*%alpha),ncol=1L)
      alpha <- c(0,1)
    }
    h <- if(type=="initial") rep(1,source$n) else
      RoCE:::calculate_glm_gradient_cpp(outcome_design,alpha,1L,1L)[,1L]
    target_h <- if(type=="initial") rep(1,target$n) else
      RoCE:::calculate_glm_gradient_cpp(target_design,alpha,1L,1L)[,1L]
    moment <- colMeans(cbind(1,target$Z_site)*target_h)
    maximum <- if(type=="initial") RoCE:::compute_lambda_max_initial_dr(source$Z_site,source$A,moment,case$arm) else
      RoCE:::compute_lambda_max_refined_dr(source$Z_site,source$A,moment,alpha,A_val=case$arm,
        W_outcome=outcome_design,calibrated=type=="calibrated",M_tau=12,truncate_initial_outcome=FALSE)
    grid <- RoCE:::build_lambda_grid(maximum,RoCE:::LAMBDA_MIN_RATIO_LOW_DIM,100L)
    set.seed(case$seed+300L)
    fold_ids <- sample(rep(1:5,length.out=sum(source$A==case$arm)))
    common <- list(Z_site=source$Z_site,A_source=source$A,lambda_grid=grid,
      n_folds=5L,max_iter=10000L,tol=1e-10,A_val=case$arm,cv_fold_id=fold_ids)
    if(type=="initial") {
      call_arguments <- c(common,list(mean_phi=moment,M_tau=12))
      function_name <- "select_lambda_cv_initial_density_ratio_cpp"
    } else {
      call_arguments <- c(common,list(mean_grad_psi=moment,alpha_init=alpha,family_int=1L,link_int=1L))
      function_name <- "select_lambda_cv_density_ratio_cpp"
      if(type=="calibrated") {
        function_name <- "select_lambda_cv_calibrated_density_ratio_cpp"
        call_arguments <- c(call_arguments,list(M_tau=12,W_outcome=outcome_design,truncate_initial_outcome=FALSE))
      }
    }
    for(solver in c("proximal_newton","coordinate_descent")) {
      name <- paste(case$config,case$p,case$arm,type,solver,sep="_")
      specifications[[name]] <- list(p=case$p,config=case$config,arm=case$arm,type=type,
        solver=solver,function_name=function_name,arguments=call_arguments,
        outcome_design=outcome_design,h=h,moment=moment)
    }
  }
}
saveRDS(specifications,file.path(output,"specifications.rds"))
script <- file.path(workflow,"validate_density_cv_certificate.R")
status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script),shQuote(check),
  shQuote(reference_library),shQuote(workflow),shQuote(output),"reference"),
  stdout=file.path(output,"reference.log"),stderr=file.path(output,"reference.err"))
stopifnot(status==0L)
}
references <- readRDS(file.path(output,"reference_results.rds"))
Sys.setenv(PKG_CPPFLAGS=paste("-I",shQuote(file.path(check,"source/src"))),PKG_CXXFLAGS="-g0 -O2")
Rcpp::sourceCpp(file.path(workflow,"profile_density_grid.cpp"),cacheDir=file.path(output,"compiled"))
rows <- list(); results <- list()
old_rows <- if(resume_mode && file.exists(file.path(output,"comparison.csv")))
  read.csv(file.path(output,"comparison.csv"),stringsAsFactors=FALSE) else NULL
check_result <- function(name,record) {
  off <- record$reference;on <- record$certificate
  stopifnot(isTRUE(all.equal(off,references[[name]]$cv,tolerance=1e-12)),
    identical(off$lambda_min,on$lambda_min),identical(off$lambda_1se,on$lambda_1se),
    isTRUE(all.equal(off,on[names(off)],tolerance=0)),
    !nzchar(record$profile_reference$error),!nzchar(record$profile_certificate$error),
    identical(is.finite(record$profile_reference$fold_scores),is.finite(record$profile_certificate$fold_scores)),
    isTRUE(all.equal(record$profile_reference$fold_scores,record$profile_certificate$fold_scores,tolerance=0)),
    isTRUE(all.equal(record$profile_reference$summary$cv_scores,off$cv_scores,tolerance=1e-12)),
    sum(record$profile_certificate$profile$certified_nonconvergence)==on$certified_nonconvergent_fold_fits,
    max(abs(as.numeric(record$fit)-as.numeric(references[[name]]$fit)))<1e-9,
    isTRUE(attr(record$fit,"converged")))
  invisible(TRUE)
}
for(name in names(specifications)) {
  specification <- specifications[[name]];a <- specification$arguments
  saved <- file.path(output,paste0(name,".rds"))
  if(resume_mode && file.exists(saved)) {
    record <- readRDS(saved)
    check_result(name,record)
    if(is.null(old_rows)) stop("Completed cases have no timing record")
    row <- old_rows[old_rows$case==name,,drop=FALSE]
    if(nrow(row)!=1L) stop("Completed case lacks an unambiguous timing record: ",name)
    rows[[name]] <- row
    cat("REUSED",name,"\n")
    next
  }
  Sys.setenv(ROCE_NUISANCE_SOLVER=specification$solver)
  RoCE:::set_nuisance_solver_cpp(specification$solver)
  native <- getFromNamespace(specification$function_name,"RoCE")
  off_time <- system.time(off <- do.call(native,c(a,list(use_kkt_certificate=FALSE))))[["elapsed"]]
  on_time <- system.time(on <- do.call(native,c(a,list(use_kkt_certificate=TRUE))))[["elapsed"]]
  profile_arguments <- list(Z=a$Z_site,A=a$A_source,linear=specification$moment,
    h_all=specification$h,lambdas=a$lambda_grid,arm=a$A_val,n_folds=5L,
    radius=if(specification$type=="refined") Inf else 12,tolerance=1e-10,max_iter=10000L,
    cv_fold_id=a$cv_fold_id)
  profile_off <- do.call(profile_density_grid_cpp,c(profile_arguments,list(use_kkt_certificate=FALSE)))
  profile_on <- do.call(profile_density_grid_cpp,c(profile_arguments,list(use_kkt_certificate=TRUE)))
  fitted <- refit(specification,on$lambda_min)
  difference <- max(abs(as.numeric(fitted)-as.numeric(references[[name]]$fit)))
  rows[[name]] <- data.frame(case=name,p=specification$p,config=specification$config,
    arm=specification$arm,type=specification$type,solver=specification$solver,
    selected_lambda=on$lambda_min,certified_fits=on$certified_nonconvergent_fold_fits,
    reference_seconds=off_time,certificate_seconds=on_time,coefficient_difference=difference)
  results[[name]] <- list(reference=off,certificate=on,profile_reference=profile_off,
    profile_certificate=profile_on,fit=fitted)
  check_result(name,results[[name]])
  temporary <- paste0(saved,".tmp")
  saveRDS(results[[name]],temporary)
  stopifnot(file.rename(temporary,saved))
  temporary <- file.path(output,"comparison.csv.tmp")
  write.csv(do.call(rbind,rows),temporary,row.names=FALSE)
  stopifnot(file.rename(temporary,file.path(output,"comparison.csv")))
  cat("PASSED",name,"certified",on$certified_nonconvergent_fold_fits,"\n")
}
stopifnot(sum(vapply(rows,function(row) row$certified_fits,numeric(1L)))>0)
writeLines(c("DENSITY_CV_PATH_AND_REFIT_PARITY_PASSED",paste("Cases:",length(rows)),
  paste("Validated:",names(rows))),
  file.path(output,"checks.txt"))

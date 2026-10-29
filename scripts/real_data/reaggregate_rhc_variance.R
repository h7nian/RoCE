#!/usr/bin/env Rscript
# Reaggregate a saved RHC analysis checkpoint; no nuisance-model fitting.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript reaggregate_rhc_variance.R ANALYSIS_RDS METHODS_CSV NEW_OUTPUT_DIR")
lib <- Sys.getenv("ROCE_DIAGNOSTIC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
stopifnot(packageVersion("RoCE") >= "0.1.1")
checkpoint <- normalizePath(args[1], mustWork = TRUE)
methods_path <- normalizePath(args[2], mustWork = TRUE)
out <- args[3]
if (file.exists(out)) stop("Use a new output directory")
dir.create(out, recursive = TRUE)
checks <- list()
check <- function(name, passed) {
  checks[[length(checks)+1L]] <<- data.frame(check = name, passed = isTRUE(passed))
  if (!isTRUE(passed)) stop("Check failed: ", name)
}
tryCatch({
  saved <- readRDS(checkpoint)
  analysis <- if (!is.null(saved$value)) saved$value else saved
  fit <- if (!is.null(analysis$tate_fit)) analysis$tate_fit else analysis
  methods <- read.csv(methods_path, stringsAsFactors = FALSE)
  original <- methods[methods$method == "RoCE", , drop = FALSE]
  check("reference_matches_saved_fit", nrow(original)==1L &&
    abs(original$estimate-fit$estimate)<1e-10 && abs(original$se-fit$se)<1e-10)
  check("joint_tate", identical(fit$aggregation_mode,"joint_tate"))
  sizes <- fit$intermediates$sample_sizes
  sources <- colnames(fit$fold_weights_by_arm$mu1)
  counts <- c(t=as.numeric(sizes$n_t),setNames(as.numeric(sizes$n_source),sources))
  check("RHC_reference_cohort",sum(counts)==5039 && counts[1L]==1698 && length(sources)==3L)
  metadata <- lapply(counts,function(n)list(n=n))
  frozen <- digest::digest(fit$arm_results,algo="sha256")
  trace("run_tate_crossfit",quote(stop("Unexpected nuisance refit")),where=asNamespace("RoCE"),print=FALSE)
  trace("run_crossfit",quote(stop("Unexpected nuisance refit")),where=asNamespace("RoCE"),print=FALSE)
  primary <- NULL; rows <- list()
  for (lambda in sort(unique(c(fit$aggregation_lambda_selection,0,.1,.25,.5,1,2,4,8,16,64,1e6)))) {
    x <- calculate_tate_crossfit_aggregation(metadata,fit$arm_results$mu1,fit$arm_results$mu0,
      lambda_selection=lambda,aggregation_mode="joint_tate",
      screening_rule=fit$aggregation_screening_rule,verbose=FALSE)
    check(paste(lambda,"leading_convention"),identical(x$variance_method,"leading") &&
      abs(x$se-x$se_fixed_weights)<1e-12)
    check(paste(lambda,"target_unchanged"),all(vapply(c("estimate","se","variance"),
      function(k)isTRUE(all.equal(x$target_only[[k]],fit$target_only[[k]],tolerance=1e-12)),logical(1))))
    if (all(x$fold_weights==0)) check(paste(lambda,"zero_weights_target_limit"),
      abs(x$estimate-x$target_only$estimate)<1e-12 && abs(x$se-x$target_only$se)<1e-12)
    if (lambda==fit$aggregation_lambda_selection) {
      check("primary_point_and_weights_unchanged",abs(x$estimate-fit$estimate)<1e-12 &&
        isTRUE(all.equal(x$fold_weights,fit$fold_weights,tolerance=1e-12)))
      check("primary_SE_matches_saved_leading",abs(x$se-fit$se_fixed_weights)<1e-12)
      check("adjusted_diagnostic_matches_previous_SE",abs(x$se_weight_adjusted-fit$se)<1e-12)
      primary <- x
    }
    rows[[length(rows)+1L]] <- data.frame(lambda=lambda,cutoff=if(lambda==0)Inf else 1/lambda,
      estimate=x$estimate,se=x$se,ci_lower=x$ci_lower,ci_upper=x$ci_upper,
      se_weight_adjusted=x$se_weight_adjusted,target_estimate=x$target_only$estimate,
      target_se=x$target_only$se,se_ratio_to_target=x$se/x$target_only$se,
      above_target=x$se-x$target_only$se>1e-10*max(1,x$target_only$se),
      all_weights_zero=all(x$fold_weights==0),nonzero_weight_coordinates=sum(x$fold_weights!=0),
      nuisance_refits=0L,variance_method="leading")
  }
  check("nuisance_scores_unchanged",identical(frozen,digest::digest(fit$arm_results,algo="sha256")))
  write.csv(do.call(rbind,rows),file.path(out,"rhc_lambda_path.csv"),row.names=FALSE)
  methods$variance_method <- "unchanged_reference"
  methods$se_weight_adjusted <- NA_real_
  i <- which(methods$method=="RoCE")
  methods$se[i] <- primary$se
  methods$ci_lower[i] <- primary$ci_lower;methods$ci_upper[i] <- primary$ci_upper
  methods$variance_method[i] <- "leading";methods$se_weight_adjusted[i] <- primary$se_weight_adjusted
  write.csv(methods,file.path(out,"rhc_methods_leading.csv"),row.names=FALSE)
  write.csv(data.frame(estimate=primary$estimate,old_se=fit$se,new_se=primary$se,
    old_ci_lower=fit$estimate-qnorm(.975)*fit$se,old_ci_upper=fit$estimate+qnorm(.975)*fit$se,
    new_ci_lower=primary$ci_lower,new_ci_upper=primary$ci_upper,
    target_se=primary$target_only$se,relative_efficiency=primary$target_only$se^2/primary$se^2,
    nuisance_refits=0L),file.path(out,"rhc_variance_change.csv"),row.names=FALSE)
  write.csv(data.frame(field=c("checkpoint_sha256","methods_sha256","package_version","data_sha256","fold_sha256"),
    value=c(digest::digest(checkpoint,algo="sha256",file=TRUE),digest::digest(methods_path,algo="sha256",file=TRUE),
      as.character(packageVersion("RoCE")),analysis$metadata$data_sha256,analysis$metadata$fold_sha256)),
    file.path(out,"provenance.csv"),row.names=FALSE)
  writeLines("SUCCESS: saved RHC fit reaggregated; zero nuisance refits.",file.path(out,"STATUS.txt"))
},error=function(e){
  writeLines(c("FAILED",conditionMessage(e)),file.path(out,"STATUS.txt"));stop(e)
},finally={
  if(length(checks))write.csv(do.call(rbind,checks),file.path(out,"verification.csv"),row.names=FALSE)
})

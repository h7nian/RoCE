#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: diagnose_n100_coverage_mechanism.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- args[[2L]]
  if (file.exists(output)) stop("diagnostic output already exists")
  source("scripts/slurm/atomic_output.R")
  source("scripts/slurm/result_provenance.R")
  summary_root <- file.path(root, "summaries/n100")
  summary_hash <- roce_sha256_file(file.path(summary_root, "sha256.txt"))
  stopifnot(summary_hash == "bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe")
  refs <- read.csv(file.path(summary_root, "input_references.csv"))
  stopifnot(identical(refs$sim_id, 10001:10100))
  all_seeds <- all_folds <- vector("list", 100L)
  for (i in 1:100) {
    directory <- file.path(root, sprintf("seed_%06d", 10000L+i))
    stopifnot(roce_sha256_file(file.path(directory,"sha256.txt")) ==
                refs$bundle_checksum_manifest_hash[i])
    manifest <- readLines(file.path(directory, "sha256.txt"))
    for (name in c("artifacts.rds", "results.csv")) {
      expected <- substr(manifest[substring(manifest,67L) == name],1L,64L)
      stopifnot(length(expected)==1L, roce_sha256_file(file.path(directory,name))==expected)
    }
    bundle <- readRDS(file.path(directory, "artifacts.rds"))
    raw <- read.csv(file.path(directory, "results.csv"))
    seed_rows <- fold_rows <- list()
    for (rho in c(0,.5,1,1.5,2,2.5)) {
      f <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
      soft <- raw[raw$rho==rho & raw$method=="one_round_crossfit_ate",]
      target <- raw[raw$rho==rho & raw$method=="target_only_ate",]
      stopifnot(nrow(soft)==1L,nrow(target)==1L)
      sizes <- f$intermediates$sample_sizes
      contributions <- numeric(2L)
      baseline <- 0
      for (k in seq_along(f$intermediates$fold_info)) {
        info <- f$intermediates$fold_info[[k]]
        raw_target <- info$varphi_ot + info$fold_target_estimate
        baseline <- baseline + sum(raw_target)/sizes$n_t
        for (j in 1:2) {
          target_part <- sum(info$zeta_components[[j]] + info$mu_pred_ts[j] - raw_target)/sizes$n_t
          source_part <- sum(info$xi_components[[j]] + info$delta_ts[j])/sizes$n_source[j]
          eta <- f$fold_weights[k,j]
          increment <- eta*(target_part+source_part)
          contributions[j] <- contributions[j]+increment
          fold_rows[[length(fold_rows)+1L]] <- data.frame(
            sim_id=10000L+i,rho,fold=k,source=paste0("s",j),weight=eta,
            wald=f$fold_wald_statistics[k,j],
            target_part_unweighted=target_part,source_part_unweighted=source_part,
            global_increment=increment)
        }
      }
      identity_error <- max(abs(baseline-target$estimate),
        abs(baseline+sum(contributions)-f$estimate), abs(f$estimate-soft$estimate))
      stopifnot(is.finite(identity_error),identity_error<1e-10)
      seed_rows[[length(seed_rows)+1L]] <- data.frame(
        sim_id=10000L+i,rho,truth=soft$truth,target=baseline,estimate=f$estimate,
        increment_s1=contributions[1],increment_s2=contributions[2],
        analytic_se=soft$se,relearn_se=soft$se_weight_relearn_bootstrap,
        identity_error=identity_error)
    }
    all_seeds[[i]] <- do.call(rbind,seed_rows)
    all_folds[[i]] <- do.call(rbind,fold_rows)
    rm(bundle,f)
    if (i%%10L==0L) cat("decomposed seeds:",i,"\n")
  }
  seeds <- do.call(rbind,all_seeds)
  folds <- do.call(rbind,all_folds)
  stopifnot(nrow(seeds)==600L,nrow(folds)==6000L)
  mcse <- function(x) sd(x)/sqrt(length(x))
  q <- qnorm(.975)
  summaries <- lapply(split(seeds,seeds$rho),function(x) {
    e <- x$estimate-x$truth
    t <- x$target-x$truth
    d <- x$increment_s1+x$increment_s2
    centered <- e-mean(e)
    empirical_sd <- sd(e)
    variance_identity_error <- abs(var(e)-var(t)-var(d)-2*cov(t,d))
    stopifnot(variance_identity_error<1e-12,
              abs(mean(e)-mean(t)-mean(x$increment_s1)-mean(x$increment_s2))<1e-12)
    weights <- folds[folds$rho==x$rho[1],]
    data.frame(rho=x$rho[1],n=nrow(x),bias=mean(e),bias_mcse=mcse(e),
      target_bias=mean(t),mean_increment_s1=mean(x$increment_s1),
      increment_s1_mcse=mcse(x$increment_s1),mean_increment_s2=mean(x$increment_s2),
      increment_s2_mcse=mcse(x$increment_s2),
      mean_weight_s1=mean(weights$weight[weights$source=="s1"]),
      mean_weight_s2=mean(weights$weight[weights$source=="s2"]),
      empirical_variance=var(e),target_variance=var(t),borrowing_variance=var(d),
      twice_target_borrowing_covariance=2*cov(t,d),variance_identity_error,
      mean_analytic_variance=mean(x$analytic_se^2),mean_relearn_variance=mean(x$relearn_se^2),
      analytic_coverage=mean(abs(e)<=q*x$analytic_se),
      analytic_misses_estimate_too_low=sum(e < -q*x$analytic_se),
      analytic_misses_estimate_too_high=sum(e > q*x$analytic_se),
      relearn_coverage=mean(abs(e)<=q*x$relearn_se),
      relearn_misses_estimate_too_low=sum(e < -q*x$relearn_se),
      relearn_misses_estimate_too_high=sum(e > q*x$relearn_se),
      oracle_centered_with_analytic_se_fraction=mean(abs(centered)<=q*x$analytic_se),
      oracle_centered_with_relearn_se_fraction=mean(abs(centered)<=q*x$relearn_se),
      oracle_empirical_sd_uncentered_fraction=mean(abs(e)<=q*empirical_sd),
      oracle_empirical_sd_centered_fraction=mean(abs(centered)<=q*empirical_sd),
      error_skewness=mean(centered^3)/mean(centered^2)^1.5,
      error_excess_kurtosis=mean(centered^4)/mean(centered^2)^2-3,
      analytic_z_mean=mean(e/x$analytic_se),analytic_z_sd=sd(e/x$analytic_se),
      relearn_z_mean=mean(e/x$relearn_se),relearn_z_sd=sd(e/x$relearn_se))
  })
  summary <- do.call(rbind,summaries)
  roce_write_atomic_directory(output,function(stage) {
    write.csv(seeds,file.path(stage,"seed_decomposition.csv"),row.names=FALSE)
    write.csv(folds,file.path(stage,"fold_source_contributions.csv"),row.names=FALSE)
    write.csv(summary,file.path(stage,"mechanism_summary.csv"),row.names=FALSE)
    writeLines(c("saved_n100_mechanism_diagnostic=complete",
      "purpose=post_checkpoint_algebra_and_oracle_thought_experiments_not_validated_CIs",
      paste0("summary_checksum_manifest_fingerprint=",summary_hash),
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/diagnose_n100_coverage_mechanism.R")),
      paste0("max_estimator_identity_error=",max(seeds$identity_error)),
      "new_mc_replications=0","primary_estimator_se_ci_changed=FALSE",
      "inference_validated=FALSE"),file.path(stage,"metadata.txt"))
    files <- list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),
               file.path(stage,"sha256.txt"))
  },caller="n100 coverage mechanism diagnostic")
  print(summary[,c("rho","bias","target_bias","mean_increment_s1","mean_increment_s2",
                    "mean_weight_s1","mean_weight_s2","analytic_coverage","relearn_coverage")],row.names=FALSE)
}

if(sys.nframe()==0L) main()

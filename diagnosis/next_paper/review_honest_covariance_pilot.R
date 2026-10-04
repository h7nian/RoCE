#!/usr/bin/env Rscript
# Reuse every prespecified fitted model; separate covariance and final means.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: repository_root completed_study new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
study <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
manifest <- read.csv(file.path(study,"manifest.csv"),stringsAsFactors=FALSE)
if(anyDuplicated(manifest$task_id) || anyDuplicated(manifest[c("scenario","K","rho","seed")])) stop("Duplicated diagnostic cases")
first <- readRDS(file.path(study,"cases",manifest$task_id[1L],"configuration.rds"))
.libPaths(c(first$library,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=first$library))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(first$library,"RoCE"))))
files <- c("candidate_score_summary.R","honest_candidate_scores.R","honest_covariance_pilot.R",
  "dispersion_bias_intervals.R","joint_dispersion_intervals.R")
for(file in files) source(file.path(repository,"diagnosis/next_paper",file))
dir.create(output,recursive=TRUE)
rows <- list(); saved_packets <- list(); parity_checks <- list()
for(index in seq_len(nrow(manifest))) {
  task <- manifest[index,]
  directory <- file.path(study,"cases",task$task_id)
  complete <- readRDS(file.path(directory,"COMPLETE.rds"))
  configuration <- readRDS(file.path(directory,"configuration.rds"))
  stopifnot(identical(complete$configuration_hash,digest::digest(configuration,algo="sha256")),
    identical(unname(complete$hashes),unname(tools::md5sum(file.path(directory,names(complete$hashes))))),
    identical(configuration$source_hashes,tools::md5sum(names(configuration$source_hashes))),
    identical(configuration$library,first$library))
  for(field in c("scenario","K","rho","seed","invalid_sources")) {
    if(configuration[[field]]!=task[[field]]) stop("Manifest differs from configuration: ",field)
  }
  set.seed(configuration$seed)
  dgp <- generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,configuration$K),p=4L,
    config=configuration$scenario,ate_deviation=configuration$rho,
    n_deviated_sites=configuration$invalid_sources,deviation_mechanism="treated_arm")
  data <- split_data_by_site(dgp)
  saved <- readRDS(file.path(directory,"fitted.rds"))
  previous <- readRDS(file.path(directory,"population_packets.rds"))
  stopifnot(identical(saved$data_hash,digest::digest(data,algo="sha256")),
    identical(saved$data_hash,previous$data_hash),
    identical(saved$configuration_hash,complete$configuration_hash),
    identical(previous$configuration_hash,complete$configuration_hash))
  fitted <- saved$fitted
  split <- split_honest_evaluation(fitted)
  variants <- list(translated_source=fitted,common_outcome=fitted)
  sources <- setdiff(fitted$sites,"t")
  for(arm in c("mu1","mu0")) for(site in sources) {
    variants$common_outcome$models[[arm]]$source[[site]]$outcome <- fitted$models[[arm]]$common$coefficients
  }
  case_packets <- list()
  for(construction in names(variants)) {
    model <- variants[[construction]]
    old <- previous$variants[[construction]]
    default <- evaluate_honest_candidates(model,data)
    stopifnot(identical(default,old$packet))
    pilot <- evaluate_honest_candidates(model,data,evaluation_indices=split$covariance)
    evaluation <- evaluate_honest_candidates(model,data,evaluation_indices=split$evaluation)
    stopifnot(all(pilot$evaluation_sizes==250),all(evaluation$evaluation_sizes==250),
      all(evaluation$training_sizes==500),identical(pilot$training_shifts,evaluation$training_shifts))
    estimated <- rescale_honest_covariance(pilot,evaluation$evaluation_sizes)
    columns <- c(1L,rep(3L,length(sources)),2L,rep(4L,length(sources)))
    population_components <- list(evaluation_sizes=old$packet$evaluation_sizes,
      target_covariance=old$population$target_covariance[columns,columns],
      common_covariance=old$population$target_covariance,
      source_private_covariances=lapply(old$population$source_moments,`[[`,"covariance"))
    population <- rescale_honest_covariance(population_components,evaluation$evaluation_sizes)
    stopifnot(max(abs(population$covariance-2*old$population$covariance))<1e-14)
    width <- nrow(evaluation$means)
    design <- cbind(mu1=rep(c(1,0),each=width),mu0=rep(c(0,1),each=width))
    selected <- which(as.numeric(old$valid)==1)
    truth <- old$population$truth[1L]-old$population$truth[2L]
    covariances <- list(conditional_population=population$covariance,
      independent_pilot=estimated$covariance,same_evaluation=evaluation$covariance)
    for(type in names(covariances)) {
      covariance <- covariances[[type]]
      oracle_valid <- joint_mean_gls(covariance,design,selected)
      coefficients <- oracle_valid$coefficients
      bias <- sum(coefficients*as.numeric(old$population$means))-truth
      variance <- drop(crossprod(coefficients,population$covariance%*%coefficients))
      inverse_cholesky <- backsolve(chol(covariance),diag(nrow(covariance)))
      relative <- crossprod(inverse_cholesky,population$covariance%*%inverse_cholesky)
      eigenvalues <- eigen(relative,symmetric=TRUE,only.values=TRUE)$values
      epsilon <- max(abs(eigenvalues-1))
      rows[[length(rows)+1L]] <- cbind(task,data.frame(construction=construction,covariance_type=type,
        evaluation_size=250L,covariance_size=if(type=="conditional_population") NA_integer_ else 250L,
        estimate=sum(coefficients*as.numeric(evaluation$means)),truth=as.numeric(truth),
        fixed_weight_bias=as.numeric(bias),fixed_weight_sd=sqrt(variance),
        bias_over_fixed_weight_sd=as.numeric(bias/sqrt(variance)),
        reported_to_fixed_weight_variance=oracle_valid$variance/variance,
        loewner_epsilon=epsilon,relative_eigenvalue_min=min(eigenvalues),relative_eigenvalue_max=max(eigenvalues),
        moment_scope=if(type=="same_evaluation") "fixed_realized_weights_only" else "conditional_on_training_and_pilot"))
    }
    case_packets[[construction]] <- list(pilot=pilot,evaluation=evaluation,
      covariance=estimated,population_covariance=population,valid=old$valid)
    parity_checks[[length(parity_checks)+1L]] <- cbind(task,data.frame(construction=construction,
      original_packet_identical=TRUE,data_hash_matches=TRUE,population_scaling_error=
        max(abs(population$covariance-2*old$population$covariance))))
  }
  saved_packets[[as.character(task$task_id)]] <- list(configuration=configuration,split=split,variants=case_packets)
  cat("Checked independent covariance split",index,"/",nrow(manifest),"\n")
}
result <- do.call(rbind,rows)
cells <- split(result,with(result,interaction(scenario,K,rho,construction,covariance_type,drop=TRUE)))
summary <- do.call(rbind,lapply(cells,function(cell) data.frame(
  scenario=cell$scenario[1L],K=cell$K[1L],rho=cell$rho[1L],construction=cell$construction[1L],
  covariance_type=cell$covariance_type[1L],training_repeats=nrow(cell),
  minimum_variance_ratio=min(cell$reported_to_fixed_weight_variance),
  median_variance_ratio=median(cell$reported_to_fixed_weight_variance),
  maximum_variance_ratio=max(cell$reported_to_fixed_weight_variance),
  median_loewner_epsilon=median(cell$loewner_epsilon),maximum_loewner_epsilon=max(cell$loewner_epsilon),
  max_absolute_standardized_bias=max(abs(cell$bias_over_fixed_weight_sd)))))
write.csv(result,file.path(output,"aggregate_diagnostics.csv"),row.names=FALSE)
write.csv(summary,file.path(output,"cell_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,parity_checks),file.path(output,"parity_checks.csv"),row.names=FALSE)
saveRDS(saved_packets,file.path(output,"pilot_packets.rds"))
saveRDS(list(study=study,manifest=manifest,manifest_md5=tools::md5sum(file.path(study,"manifest.csv")),
  source_hashes=tools::md5sum(file.path(repository,"diagnosis/next_paper",c(files,"review_honest_covariance_pilot.R")))),
  file.path(output,"provenance.rds"))
writeLines(c("HONEST_COVARIANCE_PILOT_REVIEW_PASSED",paste("Fitted repeats:",nrow(manifest)),
  "1000/site = 500 training + 250 covariance pilot + 250 final mean evaluation. No nuisance refits.",
  "The original packets are exactly reproduced; the new split is index-only and data-independent.",
  "Validity labels and population covariance remain oracle diagnostics. Five training draws/cell do not validate coverage.",
  "Same-evaluation covariance rows hold weights fixed; independent-pilot rows have conditional mean/variance interpretations.",
  "Loewner errors are realized diagnostics, not simultaneous high-probability covariance bounds."),file.path(output,"checks.txt"))

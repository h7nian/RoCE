#!/usr/bin/env Rscript
# Review every prespecified fitted diagnostic; do not select successful seeds.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: repository_root study_root new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
study <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
source(file.path(repository,"diagnosis/next_paper/honest_bias_decomposition.R"))
manifest <- read.csv(file.path(study,"manifest.csv"),stringsAsFactors=FALSE)
identity <- with(manifest,paste(scenario,K,rho,seed,sep=":"))
if(anyDuplicated(identity) || anyDuplicated(manifest$task_id)) stop("Duplicated diagnostic repeats")
candidate_rows <- list();aggregate_rows <- list();decompositions <- list()
for(index in seq_len(nrow(manifest))) {
  task <- manifest[index,]
  directory <- file.path(study,"cases",task$task_id)
  marker <- file.path(directory,"COMPLETE.rds")
  if(!file.exists(marker)) stop("Prespecified repeat is incomplete: ",task$task_id)
  complete <- readRDS(marker)
  paths <- file.path(directory,names(complete$hashes))
  stopifnot(identical(unname(complete$hashes),unname(tools::md5sum(paths))))
  configuration <- readRDS(file.path(directory,"configuration.rds"))
  configuration_hash <- digest::digest(configuration,algo="sha256")
  stopifnot(identical(complete$configuration_hash,configuration_hash),
    identical(configuration$source_hashes,tools::md5sum(names(configuration$source_hashes))))
  for(field in c("scenario","K","rho","seed","invalid_sources")) {
    if(configuration[[field]]!=task[[field]]) stop("Manifest and fitted configuration differ: ",field)
  }
  packet <- readRDS(file.path(directory,"population_packets.rds"))
  saved_fit <- readRDS(file.path(directory,"fitted.rds"))
  stopifnot(identical(packet$configuration_hash,configuration_hash),
    identical(saved_fit$configuration_hash,configuration_hash),
    identical(packet$data_hash,saved_fit$data_hash))
  candidate <- read.csv(file.path(directory,"candidate_moments.csv"),stringsAsFactors=FALSE)
  aggregate <- read.csv(file.path(directory,"oracle_valid_aggregate.csv"),stringsAsFactors=FALSE)
  candidate_rows[[index]] <- cbind(task[rep(1L,nrow(candidate)),],candidate)
  aggregate_rows[[index]] <- cbind(task[rep(1L,nrow(aggregate)),],aggregate)
  for(construction in names(packet$variants)) {
    variant <- packet$variants[[construction]]
    population <- variant$population
    result <- decompose_honest_bias(population$means,population$covariance,population$truth,variant$valid)
    fitted_drift <- sweep(population$means,2L,population$truth,"-")
    maximum <- max(abs(fitted_drift[variant$valid]))
    target_size <- unname(variant$packet$evaluation_sizes["t"])
    reference <- aggregate[aggregate$construction==construction &
      aggregate$covariance_type=="conditional_population",]
    stopifnot(nrow(reference)==1L,abs(result$tate_shift-reference$fixed_weight_bias)<1e-12,
      abs(result$oracle_sd-reference$fixed_weight_sd)<1e-12)
    decompositions[[length(decompositions)+1L]] <- cbind(task,data.frame(
      construction=construction,mu1_shift=result$mean_shift[1L],mu0_shift=result$mean_shift[2L],
      tate_shift=result$tate_shift,oracle_sd=result$oracle_sd,
      standardized_tate_shift=result$tate_shift/result$oracle_sd,
      residual_mahalanobis=result$residual_mahalanobis,
      residual_degrees_freedom=result$residual_degrees_freedom,
      residual_over_sqrt_df=result$residual_over_sqrt_df,
      level_mahalanobis=result$level_mahalanobis,projection_identity_error=result$identity_error,
      max_valid_mean_bias=maximum,scaled_max_valid_mean_bias=sqrt(target_size)*task$K^.25*maximum))
  }
}
decomposition <- do.call(rbind,decompositions)
cells <- split(decomposition,with(decomposition,interaction(scenario,K,rho,construction,drop=TRUE)))
cell_summary <- do.call(rbind,lapply(cells,function(cell) {
  data.frame(scenario=cell$scenario[1L],K=cell$K[1L],rho=cell$rho[1L],construction=cell$construction[1L],
    training_repeats=nrow(cell),median_abs_standardized_tate_shift=median(abs(cell$standardized_tate_shift)),
    max_abs_standardized_tate_shift=max(abs(cell$standardized_tate_shift)),
    median_residual_over_sqrt_df=median(cell$residual_over_sqrt_df),
    median_scaled_valid_mean_bias=median(cell$scaled_max_valid_mean_bias),
    max_scaled_valid_mean_bias=max(cell$scaled_max_valid_mean_bias))
}))
dir.create(output,recursive=TRUE)
write.csv(do.call(rbind,candidate_rows),file.path(output,"candidate_moments.csv"),row.names=FALSE)
write.csv(do.call(rbind,aggregate_rows),file.path(output,"oracle_valid_aggregate.csv"),row.names=FALSE)
write.csv(decomposition,file.path(output,"bias_decomposition.csv"),row.names=FALSE)
write.csv(cell_summary,file.path(output,"decomposition_summary.csv"),row.names=FALSE)
saveRDS(list(study=study,manifest=manifest,manifest_md5=tools::md5sum(file.path(study,"manifest.csv")),
  complete_repeats=nrow(manifest),maximum_projection_error=max(decomposition$projection_identity_error)),
  file.path(output,"review_checks.rds"))
writeLines(c("COMPLETE_HONEST_POPULATION_REVIEW_PASSED",paste("Fitted repeats:",nrow(manifest)),
  "Oracle validity and population moments; no feasible interval or coverage claim.",
  paste("Training repeats per cell:",paste(sort(unique(cell_summary$training_repeats)),collapse=", ")),
  "These are diagnostics. Equal seed labels across K do not guarantee identical target observations.",
  "Estimated-covariance GLS diagnostics hold realized weights fixed."),file.path(output,"checks.txt"))
cat("Reviewed",nrow(manifest),"complete fitted diagnostics\n")

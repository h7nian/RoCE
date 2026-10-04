#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=1L) stop("Usage: pipeline_validation_root")
root <- normalizePath(arguments[1L],mustWork=TRUE)
cases <- c("C2_p10_two_round","C3_p200_one_round")
rows <- list()
for(case in cases) {
  load_result <- function(role) {
    directory <- file.path(root,case,role)
    stopifnot(file.exists(file.path(directory,"COMPLETE")))
    readRDS(file.path(directory,"comparison.rds"))
  }
  reference <- load_result("reference")
  for(role in c("off","on")) {
    candidate <- load_result(role)
    stopifnot(identical(reference$data_hash,candidate$data_hash))
    for(field in c("scenario","p","protocol","n_per_site","K","n_folds","nlambda","seed")) {
      stopifnot(identical(reference$configuration[[field]],candidate$configuration[[field]]))
    }
    stopifnot(isTRUE(all.equal(reference$summaries,candidate$summaries,tolerance=1e-10)),
      identical(names(reference$nuisances),names(candidate$nuisances)))
    coefficient_difference <- 0
    for(name in names(reference$nuisances)) {
      old <- reference$nuisances[[name]];new <- candidate$nuisances[[name]]
      stopifnot(identical(old$cv_seed,new$cv_seed),
        isTRUE(all.equal(old,new,tolerance=1e-10)))
      coefficient_difference <- max(coefficient_difference,abs(old$coefficients-new$coefficients))
    }
    if(role=="on" && candidate$configuration$p==200L) stopifnot(candidate$certified_fits>0L)
    rows[[length(rows)+1L]] <- data.frame(case=case,role=role,nuisance_records=length(candidate$nuisances),
      maximum_coefficient_difference=coefficient_difference,certified_skips=candidate$certified_fits)
  }
}
write.csv(do.call(rbind,rows),file.path(root,"pipeline_comparison.csv"),row.names=FALSE)
writeLines(c("FULL_PIPELINE_CERTIFICATE_PARITY_PASSED",
  "1000/site,10 folds,100 lambdas,three-layer; C2p10two-round and C3p200one-round.",
  "Frozen previous-library reference versus new off/on; all three aggregation modes, arms, weights, fitted nuisance records and CV seeds checked.",
  "Separate workers may differ in hardware, load and checkpoint reuse; elapsed times are descriptive, not a controlled speed benchmark."),file.path(root,"checks.txt"))

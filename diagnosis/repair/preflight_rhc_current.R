#!/usr/bin/env Rscript
# Data/schema review only: no effect estimates or fitted nuisance models.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: checked_R_library new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
raw <- load_rhc_raw()
excluded <- c("Medicare & Medicaid","No insurance")
cohort <- build_rhc_cohort(raw,outcome="death30",site_var="ninsclas",
  site_recode=setNames(rep(NA_character_,length(excluded)),excluded))
data <- build_rhc_data_split(cohort,K=4L,target_site="Private",seed=42L)
mapping <- attr(data,"site_mapping")
counts <- do.call(rbind,lapply(names(data),function(site) {
  value <- data[[site]]
  stopifnot(identical(value$W_outcome,value$Z_site),all(is.finite(value$W_outcome)),
    all(value$A %in% 0:1),all(value$Y %in% 0:1))
  data.frame(site=site,stratum=unname(mapping[if(site=="t") "target" else site]),
    n=value$n,treated=sum(value$A==1),control=sum(value$A==0),
    working_features=ncol(value$W_outcome),
    rank_with_intercept=qr(cbind(1,value$W_outcome))$rank)
}))
pooled <- do.call(rbind,lapply(data,`[[`,"W_outcome"))
stopifnot(nrow(raw)==5735L,nrow(cohort)==5039L,ncol(pooled)==61L,
  qr(cbind(1,pooled))$rank==62L,identical(as.integer(counts$n),c(1698L,1458L,1236L,647L)))
raw_kept <- raw[!raw$ninsclas %in% excluded,,drop=FALSE]
fields <- unique(c(RoCE:::.RHC_CONTINUOUS_VARS,RoCE:::.RHC_BINARY_VARS,
                  setdiff(RoCE:::.RHC_CATEGORICAL_VARS,"ninsclas")))
missingness <- do.call(rbind,lapply(fields,function(field) {
  value <- raw_kept[[field]]
  data.frame(variable=field,missing=sum(is.na(value) | trimws(as.character(value))==""),n=length(value))
}))
support <- list()
for(site in names(data)[-1L]) for(arm in 0:1) {
  source <- data[[site]]$Z_site[data[[site]]$A==arm,,drop=FALSE]
  ranges <- apply(source,2L,range)
  constant <- ranges[2L,]-ranges[1L,] < 1e-12
  gaps <- constant & abs(colMeans(data$t$Z_site)-ranges[1L,])>1e-10
  for(column in which(gaps)) support[[length(support)+1L]] <- data.frame(site=site,arm=arm,
    feature=colnames(source)[column],source_constant=ranges[1L,column],
    target_mean=mean(data$t$Z_site[,column]),
    scope="Empirical arm-specific moment-support gap; not by itself a population nonoverlap proof")
}
dir.create(output,recursive=TRUE)
write.csv(counts,file.path(output,"cohort_sites.csv"),row.names=FALSE)
write.csv(missingness,file.path(output,"covariate_missingness.csv"),row.names=FALSE)
gaps <- if(length(support)) do.call(rbind,support) else data.frame(site=character(),arm=integer(),
  feature=character(),source_constant=numeric(),target_mean=numeric(),scope=character())
write.csv(gaps,file.path(output,"empirical_support_gaps.csv"),row.names=FALSE)
saveRDS(list(library=library_path,data_md5=tools::md5sum(rhc_csv_path()),raw_rows=nrow(raw),
  retained_rows=nrow(cohort),source_count=length(data)-1L,total_sites=length(data),
  outcome="death30",target="Private",excluded_strata=excluded,working_features=ncol(pooled),
  pooled_rank=qr(cbind(1,pooled))$rank,site_mapping=mapping,
  rhc_interface_arguments=names(formals(run_rhc_tate_experiment)),
  estimator_arguments=names(formals(run_tate_crossfit))),file.path(output,"preflight.rds"))
writeLines(c("RHC_DATA_PREFLIGHT_PASSED",
  "No nuisance models or treatment effects were fitted.",
  "The historical insurance grouping is reproduced, not newly selected from outcomes.",
  "K=4 in the RHC loader means four total groups and three sources.",
  "Cohort imputation/standardization, covariate timing and empirical support require review before a current-method fit.",
  "The historical RHC wrapper does not expose the new two-layer/joint-TATE controls."),file.path(output,"checks.txt"))
print(counts)
cat("Empirical source-arm support flags:",nrow(gaps),"\n")

# Experimental honest candidate interface. Existing source/Hou calibration
# kernels are reused from the checked installed RoCE package. Source
# candidate_score_summary.R first for score_mean_covariance(). No interval
# validity is asserted by this fitting/statistics interface.

make_honest_fold_views <- function(data_split, train_fraction=.5, calibration_folds=5L, split_seed=20260922L) {
  RoCE:::validate_algorithm_inputs(data_split,family="binomial")
  if(length(train_fraction)!=1L || !is.numeric(train_fraction) || !is.finite(train_fraction) ||
     train_fraction<=0 || train_fraction>=1 || length(calibration_folds)!=1L || !is.finite(calibration_folds) ||
     calibration_folds<2 || calibration_folds!=floor(calibration_folds) ||
     length(split_seed)!=1L || !is.finite(split_seed) || split_seed<1 || split_seed!=floor(split_seed) ||
     split_seed>.Machine$integer.max-length(data_split)) stop("Invalid honest split controls")
  sites <- c("t",setdiff(names(data_split),"t"))
  views <- lapply(seq_along(sites),function(index) {
    data <- data_split[[sites[index]]]
    if(!is.null(data$cv_group_id) && anyDuplicated(data$cv_group_id)) stop("Honest iid evaluation does not support duplicated observation origins")
    training_size <- floor(data$n*train_fraction)
    if(training_size<2L*calibration_folds || data$n-training_size<2L) stop("Insufficient observations for the requested honest split")
    # Fold membership depends only on sample count and an explicit seed.
    # Treatment/outcome stratification would change the conditional iid argument.
    permutation <- RoCE:::with_seed(as.integer(split_seed+index-1L),sample.int(data$n))
    evaluation <- sort(permutation[seq_len(data$n-training_size)])
    training <- permutation[-seq_len(data$n-training_size)]
    labels <- rep(seq_len(calibration_folds),length.out=training_size)
    indices <- c(list(evaluation),lapply(seq_len(calibration_folds),function(fold) sort(training[labels==fold])))
    folds <- lapply(indices,function(rows) list(original_idx=as.integer(rows),n=length(rows)))
    attr(folds,".data_ref") <- data
    folds
  })
  names(views) <- sites
  views
}

honest_subset <- function(data,indices) {
  result <- list(W_outcome=data$W_outcome[indices,,drop=FALSE],Z_site=data$Z_site[indices,,drop=FALSE],
                 A=data$A[indices],Y=data$Y[indices],n=length(indices),original_idx=indices)
  if(!is.null(data$cv_group_id)) result$cv_group_id <- data$cv_group_id[indices]
  result
}

fit_target_common_projection <- function(target_folds,target_calibration,A_val,nlambda,tolerance,radius,fit_cache) {
  allowed <- target_calibration$calibration_folds
  target <- RoCE:::combine_folds(target_folds,allowed)
  log_propensity <- unlist(lapply(allowed,function(fold) {
    block <- RoCE:::materialize_fold(target_folds,fold)
    initial <- target_calibration$initial_weight[[paste0("k2_",fold)]]
    predictor <- drop(cbind(1,block$Z_site) %*% initial)
    predictor <- pmax(-radius,pmin(radius,predictor))
    -log1p(exp(-predictor))
  }),use.names=FALSE)
  # The existing refined OR solver uses exp(-Z gamma) weights. Supplying
  # log(pi) with coefficient1 gives EXACT inverse-propensity weights, not odds.
  # Numerical clipping is inactive: log(pi) lies in [-radius-log(2),0].
  if(radius+log(2)>50 || max(exp(-log_propensity))>RoCE:::WEIGHT_MAX) stop("Auxiliary IPW weights exceed the solver's numerical domain")
  coefficients <- RoCE:::.fit_nuisance_training_subset("fit_unified_outcome",list(
    W_outcome=target$W_outcome,Y=target$Y,A=target$A,A_val=A_val,
    gamma_s=c(0,1),Z_site=matrix(log_propensity,ncol=1L),lambda=NULL,
    calibrated=FALSE,M_tau=radius,family_int=1L,link_int=1L,
    warm_start=as.numeric(target_calibration$outcome),nlambda=nlambda,
    max_iter=RoCE:::MAX_ITER_DEFAULT,tol=tolerance,lambda_rule="min",cv_group_id=target$cv_group_id
  ),"t_common_projection",allowed,fit_cache)
  list(coefficients=coefficients,training_folds=allowed,training_indices=target$original_idx,
       inverse_propensity=exp(-log_propensity),weight_rule="inverse_propensity")
}

fit_honest_candidates <- function(data_split,train_fraction=.5,calibration_folds=5L,
    split_seed=20260922L,nlambda=100L,nuisance_solver="proximal_newton",tolerance=1e-10,
    radius=12,communication_mode=c("one_round","two_round"),checkpoint_dir=NULL) {
  communication_mode <- match.arg(communication_mode)
  RoCE:::.validate_score_calibration_radius(radius)
  if(length(nlambda)!=1L || !is.finite(nlambda) || nlambda<2 || nlambda!=floor(nlambda) ||
     length(tolerance)!=1L || !is.finite(tolerance) || tolerance<=0) stop("Invalid nuisance grid or tolerance")
  if(!is.null(checkpoint_dir) && !startsWith(checkpoint_dir,"/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch for persistent fits")
  previous <- RoCE:::.set_nuisance_solver(nuisance_solver)
  on.exit(RoCE:::.restore_nuisance_solver(previous),add=TRUE)
  folds <- make_honest_fold_views(data_split,train_fraction,calibration_folds,split_seed)
  sites <- names(folds)
  sources <- setdiff(sites,"t")
  if(!length(sources)) stop("At least one source is required")
  allowed <- seq.int(2L,calibration_folds+1L)
  cache <- RoCE:::.new_nuisance_cache(TRUE,checkpoint_dir=checkpoint_dir)
  on.exit(RoCE:::.release_nuisance_cache(cache),add=TRUE)
  control <- list(recipe="score_derivative",target_propensity_initialization="calibrated",target_radius=radius)
  models <- lapply(c(mu1=1L,mu0=0L),function(arm) {
    target <- RoCE:::.fit_target_calibration(folds$t,allowed,arm,"binomial",radius,nlambda,
      lambda_rule="min",fit_cache=cache,layout="compact",tol=tolerance,calibration_control=control)
    common <- fit_target_common_projection(folds$t,target,arm,nlambda,tolerance,radius,cache)
    source <- lapply(sources,function(site) RoCE:::.fit_complete_source_program(
      folds$t,folds[[site]],site,allowed,communication_mode,arm,"binomial",radius,nlambda,
      RoCE:::MAX_ITER_DEFAULT,"min",cache,layout="compact",tol=tolerance,calibration_control=control))
    names(source) <- sources
    list(target=target,common=common,source=source)
  })
  partition <- lapply(folds,function(views) list(evaluation=views[[1L]]$original_idx,
    training=unlist(lapply(views[allowed],`[[`,"original_idx"),use.names=FALSE),
    calibration=lapply(views[allowed],`[[`,"original_idx")))
  hashes <- vapply(sites,function(site)
    digest::digest(honest_subset(data_split[[site]],partition[[site]]$training),algo="sha256"),character(1L))
  list(models=models,partition=partition,training_hashes=hashes,sites=sites,
    site_sizes=vapply(data_split[sites],function(site) as.integer(site$n),integer(1L)),
    configuration=list(train_fraction=train_fraction,calibration_folds=calibration_folds,split_seed=split_seed,
      nlambda=nlambda,nuisance_solver=RoCE:::nuisance_solver_cpp(),tolerance=tolerance,
      radius=radius,communication_mode=communication_mode),
    scope="Experimental honest fitted candidates; not a validated confidence-interval procedure")
}

validate_honest_evaluation_indices <- function(fitted,evaluation_indices=NULL) {
  sites <- fitted$sites
  if(is.null(evaluation_indices)) {
    evaluation_indices <- lapply(fitted$partition[sites],`[[`,"evaluation")
  }
  if(!is.list(evaluation_indices) || is.null(names(evaluation_indices)) ||
     anyDuplicated(names(evaluation_indices)) || !setequal(names(evaluation_indices),sites)) {
    stop("evaluation_indices must be a named list containing each fitted site exactly once")
  }
  indices <- lapply(sites,function(site) {
    rows <- check_score_indices(evaluation_indices[[site]],fitted$site_sizes[[site]],site)
    if(length(rows)<2L || !all(rows %in% fitted$partition[[site]]$evaluation) ||
       any(rows %in% fitted$partition[[site]]$training)) {
      stop("Evaluation indices must contain at least two held-out observations per site: ",site)
    }
    rows
  })
  setNames(indices,sites)
}

evaluate_honest_candidates <- function(fitted,data_split,keep_scores=FALSE,evaluation_indices=NULL) {
  if(!is.logical(keep_scores) || length(keep_scores)!=1L || is.na(keep_scores)) stop("keep_scores must be TRUE or FALSE")
  sites <- fitted$sites; sources <- setdiff(sites,"t"); arms <- c("mu1","mu0")
  if(!setequal(names(data_split),sites) ||
     !identical(vapply(data_split[sites],function(site) as.integer(site$n),integer(1L)),fitted$site_sizes)) stop("Evaluation site identities or sample sizes changed")
  indices <- validate_honest_evaluation_indices(fitted,evaluation_indices)
  training <- lapply(sites,function(site) honest_subset(data_split[[site]],fitted$partition[[site]]$training))
  evaluation <- lapply(sites,function(site) honest_subset(data_split[[site]],indices[[site]]))
  names(training) <- names(evaluation) <- sites
  hashes <- vapply(training,digest::digest,character(1L),algo="sha256")
  if(!identical(hashes,fitted$training_hashes)) stop("Training data changed after fitting")
  common <- anchor <- matrix(0,evaluation$t$n,2L,dimnames=list(NULL,arms))
  shifts <- matrix(0,length(sources),2L,dimnames=list(sources,arms))
  original <- adjusted <- matrix(0,length(sources)+1L,2L,dimnames=list(c("target_anchor",sources),arms))
  raw_scores <- adjusted_scores <- setNames(lapply(sources,function(site)
    matrix(0,evaluation[[site]]$n,2L,dimnames=list(NULL,arms))),sources)
  shift_identity_error <- 0
  for(arm in arms) {
    value <- if(arm=="mu1") 1L else 0L
    models <- fitted$models[[arm]]
    common[,arm] <- RoCE:::predict_glm_cpp(evaluation$t$W_outcome,models$common$coefficients,1L,1L)
    common_training <- RoCE:::predict_glm_cpp(training$t$W_outcome,models$common$coefficients,1L,1L)
    target <- RoCE:::.evaluate_target_calibration(models$target,evaluation$t,value,"binomial",fitted$configuration$radius)
    anchor[,arm] <- target$varphi_ot+target$estimate
    original["target_anchor",arm] <- adjusted["target_anchor",arm] <- mean(anchor[,arm])
    for(site in sources) {
      model <- models$source[[site]]
      source <- evaluation[[site]]
      residual <- RoCE:::calculate_correction_term_cpp(source$Z_site,source$A,source$Y,model$weight,
        model$outcome,source$W_outcome,fitted$configuration$radius,1L,1L,value)$correction_components
      prediction <- RoCE:::predict_glm_cpp(evaluation$t$W_outcome,model$outcome,1L,1L)
      prediction_training <- RoCE:::predict_glm_cpp(training$t$W_outcome,model$outcome,1L,1L)
      shifts[site,arm] <- mean(prediction_training-common_training)
      raw_scores[[site]][,arm] <- residual
      adjusted_scores[[site]][,arm] <- residual+shifts[site,arm]
      original[site,arm] <- mean(prediction)+mean(residual)
      adjusted[site,arm] <- mean(common[,arm])+mean(adjusted_scores[[site]][,arm])
      expected <- shifts[site,arm]-mean(prediction-common[,arm])
      shift_identity_error <- max(shift_identity_error,abs(adjusted[site,arm]-original[site,arm]-expected))
    }
  }
  width <- length(sites)
  target_scores <- do.call(cbind,lapply(arms,function(arm) cbind(anchor[,arm],
    matrix(common[,arm],nrow(common),length(sources)))))
  colnames(target_scores) <- c(paste0("mu1:",rownames(adjusted)),paste0("mu0:",rownames(adjusted)))
  covariance <- target_covariance <- score_mean_covariance(target_scores)
  messages <- lapply(sources,function(site) {
    raw <- raw_scores[[site]]; shifted <- adjusted_scores[[site]]; offset <- shifts[site,]
    message <- list(n=nrow(raw),raw_sum=colSums(raw),raw_square_sum=colSums(raw^2),raw_cross_sum=sum(raw[,1]*raw[,2]))
    message$adjusted_sum <- message$raw_sum+message$n*offset
    message$adjusted_square_sum <- message$raw_square_sum+2*offset*message$raw_sum+message$n*offset^2
    message$adjusted_cross_sum <- message$raw_cross_sum+offset[1]*message$raw_sum[2]+offset[2]*message$raw_sum[1]+message$n*prod(offset)
    direct <- c(colSums(shifted),colSums(shifted^2),sum(shifted[,1]*shifted[,2]))
    from_message <- c(message$adjusted_sum,message$adjusted_square_sum,message$adjusted_cross_sum)
    if(max(abs(direct-from_message))>1e-12*max(1,abs(direct))) stop("Shifted source messages do not reproduce individual scores")
    message
  })
  names(messages) <- sources
  private <- lapply(adjusted_scores,score_mean_covariance)
  for(index in seq_along(sources)) {
    columns <- c(index+1L,width+index+1L)
    covariance[columns,columns] <- covariance[columns,columns]+private[[index]]
  }
  if(shift_identity_error>1e-12*max(1,abs(adjusted))) stop("Candidate adjustment identity failed")
  result <- list(means=adjusted,original_means=original,training_shifts=shifts,covariance=covariance,
    target_covariance=target_covariance,common_covariance=score_mean_covariance(cbind(anchor,common)),
    source_private_covariances=private,source_messages=messages,
    training_sizes=vapply(training,`[[`,numeric(1L),"n"),evaluation_sizes=vapply(evaluation,`[[`,numeric(1L),"n"),
    shift_identity_error=shift_identity_error,
    scope="Empirical honest candidate statistics; nuisance bias, covariance estimation and CI calibration remain to be validated")
  if(keep_scores) result$scores <- list(target=target_scores,source=adjusted_scores,raw_source=raw_scores,common=common,anchor=anchor)
  if(!is.null(evaluation_indices)) result$evaluation_indices <- indices
  result
}

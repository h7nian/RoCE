# Reuse checked estimation kernels to separate candidate and eta-learning effects.
# No nuisance models are fitted here. The packet contains held-out scores and
# the corresponding outer-training scores saved by layer_inference_packet().

validate_calibration_packet <- function(packet) {
  sites <- names(packet$site_sizes)
  sizes <- packet$site_sizes
  if(length(sites)<2L || sites[1L]!="t" || anyDuplicated(sites) ||
     any(!is.finite(sizes)) || any(sizes<2 | sizes!=floor(sizes)) ||
     !identical(names(packet$arms),c("mu1","mu0")) ||
     !identical(names(packet$arm_truth),c("mu1","mu0"))) stop("Invalid score packet sites, sizes or arms")
  for(arm in names(packet$arms)) {
    folds <- packet$arms[[arm]]$fold_info
    if(length(folds)!=packet$n_folds || length(packet$arms[[arm]]$inner_fold_info)!=packet$n_folds) {
      stop("Score packet fold counts differ")
    }
    for(index in seq_along(sites)) {
      ids <- unlist(lapply(folds,function(fold) if(index==1L) fold$target_idx else fold$source_idx[[index-1L]]),use.names=FALSE)
      if(any(!is.finite(ids)) || any(ids!=floor(ids)) ||
         !identical(sort(as.integer(ids)),seq_len(sizes[[index]]))) stop("Outer score indices must cover each site exactly once")
    }
  }
  invisible(packet)
}

validate_paired_score_packets <- function(first,second) {
  validate_calibration_packet(first)
  validate_calibration_packet(second)
  if(is.null(first$data_hash) || length(first$data_hash)!=1L ||
     !grepl("^[0-9a-f]{64}$",first$data_hash) || !identical(first$data_hash,second$data_hash) ||
     !identical(first$site_sizes,second$site_sizes) || !identical(first$arm_truth,second$arm_truth) ||
     !identical(first$n_folds,second$n_folds)) stop("Paired score packets must use identical data and truth")
  for(arm in names(first$arms)) for(fold in seq_len(first$n_folds)) {
    columns <- c("target_idx","source_idx")
    if(!identical(first$arms[[arm]]$fold_info[[fold]][columns],
                  second$arms[[arm]]$fold_info[[fold]][columns])) stop("Paired score packets must use identical outer folds")
  }
  invisible(TRUE)
}

fit_packet_joint_weights <- function(packet,aggregation_cutoff=2) {
  validate_calibration_packet(packet)
  if(length(aggregation_cutoff)!=1L || !is.numeric(aggregation_cutoff) ||
     !is.finite(aggregation_cutoff) || aggregation_cutoff<=0) stop("aggregation_cutoff must be positive")
  if(!identical(packet$aggregation_screening_rule,"soft_penalty")) stop("The diagnostic currently requires the soft Wald penalty")
  count <- length(packet$site_sizes)-1L
  weights <- list(mu1=matrix(0,packet$n_folds,count),mu0=matrix(0,packet$n_folds,count))
  moments <- fits <- records <- vector("list",packet$n_folds)
  for(fold in seq_len(packet$n_folds)) {
    records[[fold]] <- RoCE:::.joint_inner_records(packet$arms$mu1$inner_fold_info[[fold]],
                                                  packet$arms$mu0$inner_fold_info[[fold]])
    for(site in seq_along(packet$site_sizes)) {
      ids <- unlist(records[[fold]]$ids[[site]],use.names=FALSE)
      outer <- packet$arms$mu1$fold_info[[fold]]
      held_out <- if(site==1L) outer$target_idx else outer$source_idx[[site-1L]]
      allowed <- setdiff(seq_len(packet$site_sizes[[site]]),held_out)
      if(any(!is.finite(ids)) || any(ids!=floor(ids)) ||
         !identical(sort(as.integer(ids)),as.integer(allowed))) stop("Eta scores must cover only the outer-training observations")
    }
    moments[[fold]] <- RoCE:::.joint_inner_moments(records[[fold]])
    fits[[fold]] <- RoCE:::.solve_joint_weights(moments[[fold]],1/aggregation_cutoff,"soft_penalty")
    weights$mu1[fold,] <- fits[[fold]]$weights[seq_len(count)]
    weights$mu0[fold,] <- fits[[fold]]$weights[count+seq_len(count)]
  }
  list(packet=packet,weights=weights,moments=moments,fits=fits,records=records,
       aggregation_cutoff=aggregation_cutoff)
}

evaluate_packet_weights <- function(packet,weights,weight_learning=NULL) {
  validate_calibration_packet(packet)
  sites <- names(packet$site_sizes); sizes <- packet$site_sizes
  sources <- sites[-1L]; count <- length(sources); total <- sum(sizes)
  if(!identical(names(weights),c("mu1","mu0")) || any(vapply(weights,function(value)
      !is.matrix(value) || !is.numeric(value) || !identical(dim(value),c(as.integer(packet$n_folds),count)) ||
      any(!is.finite(value)),logical(1L)))) stop("Weights must have one row per fold and one column per source, for each arm")
  if(!is.null(weight_learning)) {
    validate_paired_score_packets(packet,weight_learning$packet)
    if(!identical(weights,weight_learning$weights)) stop("Weights differ from their saved learning program")
  }
  raw <- lapply(names(packet$arms),function(arm) RoCE:::.compute_phase3_all_phi(
    packet$n_folds,weights[[arm]],packet$arms[[arm]]$fold_info,sizes[1L],sizes[-1L],
    total,sources,count,FALSE))
  names(raw) <- names(packet$arms)
  direct <- lapply(raw,function(value) {
    groups <- split(value,rep(seq_along(sizes),sizes))
    unlist(lapply(groups,function(group) (group-mean(group))/total),use.names=FALSE)
  })
  indirect <- lapply(direct,function(value) value*0)
  if(!is.null(weight_learning)) for(arm in names(packet$arms)) {
    increments <- RoCE:::.outer_fold_increments(packet$arms[[arm]]$fold_info,sizes[1L],sizes[-1L])
    gradients <- lapply(sizes,function(size) numeric(size))
    positions <- if(arm=="mu1") seq_len(count) else count+seq_len(count)
    for(fold in seq_len(packet$n_folds)) {
      increment <- numeric(2L*count); increment[positions] <- increments[fold,]
      derivative <- RoCE:::.joint_weight_derivative(weight_learning$moments[[fold]],
        weight_learning$fits[[fold]],"soft_penalty",increment)
      for(site in seq_along(sizes)) {
        ids <- unlist(weight_learning$records[[fold]]$ids[[site]],use.names=FALSE)
        gradients[[site]][ids] <- gradients[[site]][ids]+derivative$sites[[site]]
      }
    }
    indirect[[arm]] <- unlist(gradients,use.names=FALSE)
  }
  fixed <- crossprod(do.call(cbind,direct))
  covariance <- crossprod(do.call(cbind,Map(`+`,direct,indirect)))
  means <- vapply(raw,mean,numeric(1L))
  contrast <- c(1,-1)
  list(estimates=c(means,tate=means[["mu1"]]-means[["mu0"]]),
    covariance=covariance,covariance_fixed_weights=fixed,
    tate_variance=drop(crossprod(contrast,covariance%*%contrast)),
    tate_variance_fixed_weights=drop(crossprod(contrast,fixed%*%contrast)),
    scope="Score diagnostic; eta sensitivity holds nuisance fits and active sets fixed, not a full-refit variance")
}

summarize_packet_candidates <- function(packet) {
  validate_calibration_packet(packet)
  sources <- names(packet$site_sizes)[-1L]
  candidates <- c("target_anchor",sources)
  result <- lapply(seq_along(candidates),function(index) {
    weights <- matrix(0,packet$n_folds,length(sources))
    if(index>1L) weights[,index-1L] <- 1
    value <- evaluate_packet_weights(packet,list(mu1=weights,mu0=weights))
    data.frame(candidate=candidates[index],quantity=c("mu1","mu0","tate"),
      estimate=as.numeric(value$estimates),
      truth=c(packet$arm_truth,packet$arm_truth[["mu1"]]-packet$arm_truth[["mu0"]]),
      score_variance=c(diag(value$covariance_fixed_weights),value$tate_variance_fixed_weights),
      covariance_mu1_mu0=value$covariance_fixed_weights[1L,2L])
  })
  do.call(rbind,result)
}

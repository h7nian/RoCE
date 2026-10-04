# Independent covariance-pilot support for the experimental honest candidates.
# Source candidate_score_summary.R and honest_candidate_scores.R first.
# This does not change the current paper's cross-fitted estimator.

split_honest_evaluation <- function(fitted,pilot_fraction=.5) {
  if(!is.numeric(pilot_fraction) || length(pilot_fraction)!=1L ||
     !is.finite(pilot_fraction) || pilot_fraction<=0 || pilot_fraction>=1) {
    stop("pilot_fraction must lie strictly between zero and one")
  }
  held_out <- validate_honest_evaluation_indices(fitted)
  pilot_sizes <- floor(lengths(held_out)*pilot_fraction)
  if(any(pilot_sizes<2L | lengths(held_out)-pilot_sizes<2L)) {
    stop("Each site needs at least two covariance-pilot and two evaluation observations")
  }
  # Prespecified original-index order only: no treatment/outcome balancing,
  # no extra randomness, and no use of fitted values or observed score values.
  covariance <- evaluation <- held_out
  for(site in fitted$sites) {
    rows <- sort(held_out[[site]])
    pilot <- seq_len(pilot_sizes[[site]])
    covariance[[site]] <- rows[pilot]
    evaluation[[site]] <- rows[-pilot]
  }
  list(covariance=covariance,evaluation=evaluation)
}

rescale_honest_covariance <- function(packet,evaluation_sizes) {
  pilot_sizes <- packet$evaluation_sizes
  sites <- names(pilot_sizes)
  sources <- setdiff(sites,"t")
  valid_sizes <- function(values) {
    is.numeric(values) && !is.null(names(values)) &&
      !anyDuplicated(names(values)) && setequal(names(values),sites) &&
      all(is.finite(values)) && all(values>=2 & values==floor(values))
  }
  if(!length(sources) || !identical(sites,c("t",sources)) ||
     !valid_sizes(pilot_sizes) || !valid_sizes(evaluation_sizes)) {
    stop("Covariance and evaluation sizes must name every site with integer counts of at least two")
  }
  evaluation_sizes <- evaluation_sizes[sites]
  ratios <- pilot_sizes/evaluation_sizes
  if(!is.list(packet$source_private_covariances) ||
     !identical(names(packet$source_private_covariances),sources)) {
    stop("Private covariance blocks must follow the packet's source order")
  }
  check_matrix <- function(value,size) {
    if(!is.matrix(value) || !is.numeric(value) || !identical(dim(value),c(size,size)) ||
       any(!is.finite(value)) || max(abs(value-t(value)))>1e-12*max(1,abs(value))) {
      stop("Invalid covariance block dimensions, values or symmetry")
    }
    value
  }
  width <- length(sites)
  target <- check_matrix(packet$target_covariance,2L*width)*ratios[["t"]]
  common <- check_matrix(packet$common_covariance,4L)*ratios[["t"]]
  private <- lapply(sources,function(site) {
    check_matrix(packet$source_private_covariances[[site]],2L)*ratios[[site]]
  })
  names(private) <- sources
  covariance <- target
  for(index in seq_along(sources)) {
    columns <- c(index+1L,width+index+1L)
    covariance[columns,columns] <- covariance[columns,columns]+private[[index]]
  }
  # score_mean_covariance uses centered cross-products / n_pilot^2.
  # Multiplication by n_pilot/n_eval keeps that empirical-moment convention;
  # it does not silently introduce a Bessel correction or a ridge adjustment.
  list(covariance=covariance,target_covariance=target,common_covariance=common,
    source_private_covariances=private,covariance_sizes=pilot_sizes,
    evaluation_sizes=evaluation_sizes)
}

# Exact likelihood calculation for randomized binary outcomes. No Gaussian
# approximation is needed for the lower-bound experiment.
bernoulli_mixture_chisquare <- function(n_per_site, standardized_shift, valid_probability) {
  if(length(n_per_site)!=1L || !is.finite(n_per_site) || n_per_site<1 || n_per_site!=floor(n_per_site) ||
     length(standardized_shift)!=1L || !is.finite(standardized_shift) || standardized_shift<0 ||
     length(valid_probability)!=1L || !is.finite(valid_probability) || valid_probability<=.5 || valid_probability>=1) stop("Invalid likelihood-mixture arguments")
  if(n_per_site==1L || standardized_shift==0) return(0)
  ratio <- valid_probability/(1-valid_probability)
  degree <- 2:n_per_site
  large <- log1p(-valid_probability)+degree*log(ratio)
  relative <- exp(log(valid_probability)-large)
  log_moment <- large+ifelse(degree%%2L==0L,log1p(relative),log1p(-relative))
  terms <- lchoose(n_per_site,degree)+degree*(2*log(standardized_shift)-log(n_per_site))+2*log_moment
  sum(exp(terms))
}

bernoulli_adaptation_lower_bound <- function(n_per_site, num_sources, valid_fraction,
    outcome_probability=.6, treatment_probability=.5, alpha=.05, shift_constant=.3) {
  if(length(num_sources)!=1L || !is.finite(num_sources) || num_sources<1 || num_sources!=floor(num_sources) ||
     any(!is.finite(c(valid_fraction,outcome_probability,treatment_probability,alpha,shift_constant))) ||
     any(c(valid_fraction,outcome_probability,treatment_probability,alpha)<=0) ||
     any(c(valid_fraction,outcome_probability,treatment_probability)>=1) || alpha>=.5 || shift_constant<=0) stop("Invalid lower-bound model parameters")
  probability <- (1+valid_fraction)/2
  ratio <- probability/(1-probability)
  standardized <- shift_constant/sqrt(ratio)*num_sources^(-1/4)
  per_source <- bernoulli_mixture_chisquare(n_per_site,standardized,probability)
  standard_error <- sqrt(outcome_probability*(1-outcome_probability)/(n_per_site*treatment_probability))
  shift <- standard_error*standardized
  if(outcome_probability-ratio*shift<=0 || outcome_probability+shift>=1) stop("The local alternative leaves the Bernoulli probability range")
  divergence <- expm1(n_per_site*log1p(standardized^2/n_per_site)+num_sources*log1p(per_source))
  conditioning <- pbinom(ceiling(valid_fraction*num_sources)-1L,num_sources,probability)
  total_variation <- min(1,.5*sqrt(divergence)+conditioning)
  data.frame(n_per_site=n_per_site,num_sources=num_sources,valid_fraction=valid_fraction,
    shift=shift,source_chisquare=per_source,joint_chisquare=divergence,
    conditioning_cost=conditioning,total_variation_bound=total_variation,
    expected_length_lower_bound=shift*max(0,1-2*alpha-total_variation))
}

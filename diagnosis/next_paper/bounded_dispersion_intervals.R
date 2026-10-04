# Independent bounded observations within each site. This does not assume
# Gaussian errors or known outcome variances; it does assume the supplied
# range and validity count. Source joint_dispersion_intervals.R for the
# common interval-intersection helper.

bounded_site_dispersion <- function(site_sums, site_squares, n_per_site, score_range) {
  site_sums <- as.matrix(site_sums); site_squares <- as.matrix(site_squares)
  if(!is.numeric(site_sums) || !is.numeric(site_squares) || !identical(dim(site_sums),dim(site_squares)) ||
     nrow(site_sums)<1L || ncol(site_sums)<2L ||
     any(!is.finite(c(site_sums,site_squares))) || length(n_per_site)!=1L || !is.finite(n_per_site) ||
     n_per_site<2 || n_per_site!=floor(n_per_site) || length(score_range)!=1L ||
     !is.finite(score_range) || score_range<=0) stop("Invalid bounded-site summaries")
  if(any(site_squares < site_sums^2/n_per_site-1e-10*pmax(1,site_squares))) stop("Inconsistent sums and squared sums")
  count <- ncol(site_sums)
  means <- site_sums/n_per_site
  within <- rowMeans((site_sums^2-site_squares)/(n_per_site*(n_per_site-1)))
  between <- (rowSums(means)^2-rowSums(means^2))/(count*(count-1))
  frobenius_squared <- 1/(count*n_per_site*(n_per_site-1)) + 1/(count*(count-1)*n_per_site^2)
  variance_upper <- score_range^2/4
  list(statistic=within-between,pooled_mean=rowMeans(means),
       linear_variance_factor=4*variance_upper/((count-1)*n_per_site),
       null_variance_bound=2*variance_upper^2*frobenius_squared,
       frobenius_squared=frobenius_squared,operator_norm=1/((count-1)*n_per_site),
       num_sources=count,n_per_site=n_per_site,score_range=score_range)
}

cantelli_dispersion_upper <- function(dispersion,failure_probability) {
  if(length(failure_probability)!=1L || !is.finite(failure_probability) || failure_probability<=0 || failure_probability>=1) stop("Failure probability must lie in (0,1)")
  multiplier <- (1-failure_probability)/failure_probability
  linear <- multiplier*dispersion$linear_variance_factor
  constant <- multiplier*dispersion$null_variance_bound
  observed <- dispersion$statistic
  # If the discriminant is negative the one-sided concentration event has
  # no feasible positive B. A nonnegative output is still conservative on
  # that event; on the concentration event the upper quadratic root applies.
  pmax(0,observed+linear/2+sqrt(pmax(0,linear*observed+linear^2/4+constant)))
}

bounded_source_interval <- function(target_sums,n_target,site_sums,site_squares,n_per_site,
                                    score_range,valid_minimum,alpha=.05,reference_fraction=.5,dispersion_fraction=.5) {
  dispersion <- bounded_site_dispersion(site_sums,site_squares,n_per_site,score_range)
  count <- dispersion$num_sources
  if(length(n_target)!=1L || !is.finite(n_target) || n_target<1 || n_target!=floor(n_target) ||
     length(target_sums)!=length(dispersion$pooled_mean) || any(!is.finite(target_sums)) ||
     length(valid_minimum)!=1L || !is.finite(valid_minimum) || valid_minimum<0 ||
     valid_minimum>count || valid_minimum!=floor(valid_minimum)) stop("Invalid target summaries or valid-source count")
  probabilities <- c(alpha,reference_fraction,dispersion_fraction)
  if(any(vapply(list(alpha,reference_fraction,dispersion_fraction),length,integer(1L))!=1L) ||
     any(!is.finite(probabilities)) || any(probabilities<=0 | probabilities>=1)) stop("Error allocations must lie in (0,1)")
  target <- target_sums/n_target
  if(valid_minimum==0L || valid_minimum==count) {
    total <- if(valid_minimum==0L) n_target else n_target+count*n_per_site
    estimate <- if(valid_minimum==0L) target else (target_sums+rowSums(as.matrix(site_sums)))/total
    radius <- score_range*sqrt(log(2/alpha)/(2*total))
    return(data.frame(estimate=estimate,lower=estimate-radius,upper=estimate+radius,
                     dispersion_upper=0,empty_intersection=FALSE))
  }
  source_alpha <- alpha*(1-reference_fraction)
  upper <- cantelli_dispersion_upper(dispersion,source_alpha*dispersion_fraction)
  bias_factor <- (count-valid_minimum)*(count-1)/(valid_minimum*count)
  noise_radius <- score_range*sqrt(log(2/(source_alpha*(1-dispersion_fraction)))/(2*count*n_per_site))
  source_radius <- noise_radius+sqrt(bias_factor*upper)
  target_radius <- score_range*sqrt(log(2/(alpha*reference_fraction))/(2*n_target))
  interval <- intersect_reference_interval(dispersion$pooled_mean-source_radius,dispersion$pooled_mean+source_radius,
                                          target-target_radius,target+target_radius,"shorter")
  data.frame(estimate=(interval$lower+interval$upper)/2,lower=interval$lower,upper=interval$upper,
             dispersion_upper=upper,empty_intersection=interval$empty_intersection)
}

bounded_tate_interval <- function(target_sums,n_target,site_sums,site_squares,n_per_site,
                                  score_ranges,valid_minimum,alpha=.05) {
  arms <- c("mu1","mu0")
  target_sums <- as.matrix(target_sums)
  if(!identical(colnames(target_sums),arms) || !identical(names(site_sums),arms) ||
     !identical(names(site_squares),arms) || !identical(names(score_ranges),arms) ||
     !identical(names(valid_minimum),arms) ||
     !identical(dim(as.matrix(site_sums$mu1)),dim(as.matrix(site_sums$mu0))) ||
     length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1) stop("Provide aligned mu1/mu0 summaries and named arm controls")
  intervals <- lapply(arms,function(arm)
    bounded_source_interval(target_sums[,arm],n_target,site_sums[[arm]],site_squares[[arm]],
      n_per_site,score_ranges[[arm]],valid_minimum[[arm]],alpha=alpha/2))
  names(intervals) <- arms
  lower <- intervals$mu1$lower-intervals$mu0$upper
  upper <- intervals$mu1$upper-intervals$mu0$lower
  list(estimate=(lower+upper)/2,lower=lower,upper=upper,arm_intervals=intervals,
       scope="Bounded independent site observations; cross-arm dependence is allowed by the union bound")
}

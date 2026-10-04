# Finite-population, model-based RHC sensitivity experiments. These generators
# have known causal truth; they do not estimate bias in the observed RHC cohort.
rhc_bounded_logistic <- function(coefficients, design, probabilities, radius) {
  desired_mean <- sum(probabilities * plogis(drop(design %*% coefficients)))
  stopifnot(abs(qlogis(desired_mean)) < radius)
  slopes <- drop(design[, -1L, drop = FALSE] %*% coefficients[-1L])
  candidate <- function(scale) {
    intercept <- uniroot(function(value) sum(probabilities * plogis(value + scale * slopes)) -
      desired_mean, c(-50 - max(abs(slopes)), 50 + max(abs(slopes))), tol = 1e-12)$root
    list(coefficients = c(intercept, scale * coefficients[-1L]), scale = scale,
      max_logit = max(abs(intercept + scale * slopes)))
  }
  if (candidate(1)$max_logit <= radius) return(candidate(1))
  lower <- 0; upper <- 1
  for (iteration in seq_len(55L)) {
    middle <- (lower + upper) / 2
    if (candidate(middle)$max_logit <= radius) lower <- middle else upper <- middle
  }
  candidate(lower)
}

rhc_joint_tilt <- function(design, target_probability, slopes, arm_probability, radius) {
  stopifnot(length(arm_probability) == 2L, all(arm_probability > 0),
    abs(sum(arm_probability) - 1) < 1e-12, max(-log(arm_probability)) < radius)
  linear <- design[, -1L, drop = FALSE] %*% slopes
  candidate <- function(scale) {
    tilt <- exp(scale * linear)
    normalizer <- colSums(target_probability * tilt)
    joint <- sweep(target_probability * tilt, 2L, arm_probability / normalizer, `*`)
    log_weight <- sweep(-scale * linear, 2L, log(normalizer / arm_probability), `+`)
    list(joint = joint, log_weight = log_weight, scale = scale,
      coefficients = rbind(log(arm_probability / normalizer), scale * slopes))
  }
  lower <- 0; upper <- 1
  if (max(abs(candidate(1)$log_weight)) <= radius) return(candidate(1))
  for (iteration in seq_len(55L)) {
    middle <- (lower + upper) / 2
    if (max(abs(candidate(middle)$log_weight)) <= radius) lower <- middle else upper <- middle
  }
  candidate(lower)
}

rhc_check_population <- function(population, design) {
  target <- population$sites$t
  for (site in population$sites) {
    stopifnot(length(site$probability) == nrow(design), all(is.finite(site$probability)),
      all(site$probability > 0), abs(sum(site$probability) - 1) < 1e-10,
      all(site$propensity > 0 & site$propensity < 1),
      identical(dim(site$outcome), c(nrow(design), 2L)),
      all(site$outcome > 0 & site$outcome < 1))
  }
  for (name in setdiff(names(population$sites), "t")) {
    site <- population$sites[[name]]
    joint <- site$probability * cbind(1 - site$propensity, site$propensity)
    weight <- target$probability / joint
    stopifnot(max(abs(joint * weight - target$probability)) < 1e-12)
    if (!is.null(site$weight_coefficients)) {
      stopifnot(max(abs(log(weight) + design %*% site$weight_coefficients)) < 1e-10,
        max(abs(log(weight))) <= population$weight_radius + 1e-10)
    }
  }
  truth <- sum(target$probability * (target$outcome[, 2L] - target$outcome[, 1L]))
  stopifnot(abs(truth - population$truth) < 1e-12)
  invisible(TRUE)
}

rhc_generate_population_sample <- function(template, scenario, seed) {
  stopifnot(scenario %in% names(template$populations), seed > 0L,
    seed == as.integer(seed))
  population <- template$populations[[scenario]]
  set.seed(seed)
  data <- template$observed_structure
  oracle <- list()
  for (name in names(data)) {
    site <- population$sites[[name]]
    index <- sample.int(nrow(template$features), data[[name]]$n, replace = TRUE,
      prob = site$probability)
    treatment <- rbinom(length(index), 1L, site$propensity[index])
    outcome <- rbinom(length(index), 1L, site$outcome[cbind(index, treatment + 1L)])
    data[[name]]$X <- data[[name]]$X_dagger <- template$raw_features[index, , drop = FALSE]
    data[[name]]$W_outcome <- data[[name]]$Z_site <- template$features[index, , drop = FALSE]
    data[[name]]$A <- treatment; data[[name]]$Y <- outcome
    # Repeated profiles represent independent new patients, not bootstrap copies
    # of one observed outcome. Therefore they must not be grouped in nuisance CV.
    data[[name]]$cv_group_id <- NULL
    if (name == "t") oracle <- list(index = index, outcome = site$outcome[index, , drop = FALSE],
      propensity = site$propensity[index])
  }
  m0 <- oracle$outcome[, 1L]; m1 <- oracle$outcome[, 2L]
  e <- oracle$propensity; a <- data$t$A; y <- data$t$Y
  score <- m1 - m0 + a * (y - m1) / e - (1 - a) * (y - m0) / (1 - e)
  target <- population$sites$t
  delta <- target$outcome[, 2L] - target$outcome[, 1L]
  conditional_variance <- target$outcome[, 2L] * (1 - target$outcome[, 2L]) / target$propensity +
    target$outcome[, 1L] * (1 - target$outcome[, 1L]) / (1 - target$propensity)
  variance <- sum(target$probability * ((delta - population$truth)^2 + conditional_variance)) / data$t$n
  list(data = data, truth = population$truth,
    oracle = c(estimate = mean(score), se = sqrt(var(score) / length(score)), population_se = sqrt(variance)))
}

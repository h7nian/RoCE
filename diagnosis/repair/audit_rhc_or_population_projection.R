#!/usr/bin/env Rscript
# Numerically check the population calibration limits in the clean OR branch.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
template <- readRDS(arguments[1L]); output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
population <- template$populations$O_moderate_mix
stopifnot(!is.null(population))
design <- cbind(1, template$features)
target <- population$sites$t
radius <- 3
records <- list(); fits <- list()
for (site_name in setdiff(names(population$sites), "t")) for (arm in 0:1) {
  site <- population$sites[[site_name]]
  joint <- site$probability * if (arm == 1L) site$propensity else 1 - site$propensity
  for (stage in c("initial", "outcome_derivative")) {
    h <- if (stage == "initial") rep(1, nrow(design)) else
      target$outcome[, arm + 1L] * (1 - target$outcome[, arm + 1L])
    weight <- joint * h
    moment <- colSums(design * (target$probability * h))
    objective <- function(coefficient) {
      predictor <- drop(design %*% coefficient)
      clipped <- pmax(-radius, pmin(radius, predictor))
      sum(moment * coefficient) + sum(weight * exp(-clipped) * (1 + clipped - predictor))
    }
    gradient <- function(coefficient) {
      predictor <- drop(design %*% coefficient)
      moment - drop(crossprod(design, weight * exp(-pmax(-radius, pmin(radius, predictor)))))
    }
    start <- c(log(sum(weight) / moment[1L]), rep(0, ncol(design) - 1L))
    fit <- optim(start, objective, gradient, method = "BFGS",
      control = list(maxit = 5000L, reltol = 1e-14))
    predictor <- drop(design %*% fit$par)
    q <- exp(-pmax(-radius, pmin(radius, predictor)))
    residual <- max(abs(gradient(fit$par)))
    hessian <- crossprod(design, design * (weight * q * (abs(predictor) < radius)))
    eigenvalues <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
    stopifnot(fit$convergence == 0L, residual < 2e-7, min(eigenvalues) > 1e-9)
    key <- paste(site_name, arm, stage, sep = ":")
    records[[key]] <- data.frame(site = site_name, arm = arm, stage = stage,
      gradient_max = residual, hessian_minimum_eigenvalue = min(eigenvalues),
      minimum_clipping_margin = min(abs(abs(predictor) - radius)),
      predictor_minimum = min(predictor), predictor_maximum = max(predictor),
      merged_weight_mass = sum(joint * q))
    fits[[key]] <- list(coefficients = fit$par, gradient = gradient(fit$par),
      eigenvalues = eigenvalues, optimizer = fit[c("value", "counts", "convergence")])
  }
}
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, records), file.path(output, "population_limits.csv"), row.names = FALSE)
saveRDS(fits, file.path(output, "population_coefficients.rds"))
writeLines("All12 initial/calibrated population moment systems have small residuals andpositive local Hessians at radius3; finite-sample inference still requires empirical validation.",
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, records))

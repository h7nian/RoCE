#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: check_dgp_structure.R FULL_FITS_RDS OUTPUT_PREFIX")
.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
artifact <- readRDS(args[1L])
p <- ncol(artifact$data_split$t$W_outcome) %/% 2L
n_signal <- RoCE:::FACE_SIGNAL_COORDINATES
stopifnot(p > n_signal, ncol(artifact$data_split$t$W_outcome) == 2L * p)
noise_coordinates <- seq.int(n_signal + 1L, p)
noise_coefficients <- c(noise_coordinates, p + noise_coordinates)
support_rows <- list()
for (arm_name in c("mu1", "mu0")) {
  arm <- artifact$reference$arm_results[[arm_name]]
  for (k1 in seq_len(arm$n_folds)) for (site in names(arm$fold_results[[k1]]$source_results)) {
    fit <- arm$fold_results[[k1]]$source_results[[site]]
    for (model in c("alpha_ts", "gamma_s")) {
      slopes <- as.numeric(fit[[model]])[-1L]
      support_rows[[length(support_rows) + 1L]] <- data.frame(
        arm = arm_name, k1, site, model,
        active_coefficients = sum(abs(slopes) > 1e-6),
        active_outcome_irrelevant_coefficients = sum(abs(slopes[noise_coefficients]) > 1e-6)
      )
    }
  }
}
support <- do.call(rbind, support_rows)
write.csv(support, paste0(args[2L], "_support.csv"), row.names = FALSE)

# For independent SN(kappa,1,nu) coordinates and normal target coordinates,
# log(f_source/f_target) = sum_l log(2 * Phi(nu * (X_l-kappa))). Every coordinate
# enters for nu>0, including the 96 coordinates absent from the outcome/PS signals.
centered_x <- artifact$data_split$t$W_outcome[, seq_len(p), drop = FALSE]
source_count <- length(artifact$folds$source_folds)
skewness <- seq(RoCE:::FACE_NU_SOURCE_MAX / source_count,
                RoCE:::FACE_NU_SOURCE_MAX, length.out = source_count)
density <- do.call(rbind, lapply(seq_len(source_count), function(source_index) {
  terms <- log(2) + pnorm(skewness[source_index] * centered_x, log.p = TRUE)
  data.frame(source_index, skewness = skewness[source_index], p,
             outcome_signal_coordinates = n_signal,
             shifted_coordinates = p,
             sd_log_density_ratio = sd(rowSums(terms)),
             sd_log_density_ratio_from_noise_coordinates =
               sd(rowSums(terms[, noise_coordinates, drop = FALSE])))
}))
write.csv(density, paste0(args[2L], "_density.csv"), row.names = FALSE)
print(aggregate(cbind(active_coefficients, active_outcome_irrelevant_coefficients) ~
                  model + arm, support, mean), row.names = FALSE)
print(density, row.names = FALSE)

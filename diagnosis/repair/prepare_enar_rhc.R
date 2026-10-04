#!/usr/bin/env Rscript
# Export the explicitly selected Private/min application without refitting models.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
library_path <- normalizePath(arguments[1L], mustWork = TRUE)
selection <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
provenance <- jsonlite::fromJSON(file.path(selection, "provenance.json"))
for (path in names(provenance)) stopifnot(identical(
  digest::digest(path, file = TRUE, algo = "sha256"), provenance[[path]]))
configuration <- jsonlite::fromJSON(file.path(selection, "configuration.json"))
if (!is.null(configuration$collation)) {
  stopifnot(identical(Sys.setlocale("LC_COLLATE", configuration$collation), configuration$collation))
}
stopifnot(configuration$report_target_site == "Private", configuration$report_nuisance_lambda_rule == "min",
          is.null(configuration$source_lambda_rules), configuration$crossfit_layers == 2L,
          configuration$aggregation_mode == "joint_tate", configuration$aggregation_cutoff == 2)
primary <- configuration$primary_fit
fit <- readRDS(file.path(primary, "checkpoints/analysis.rds"))$value
metadata <- fit$metadata
preprocessing <- metadata$preprocessing
if (is.null(preprocessing)) preprocessing <- "cohort"
if (!is.null(configuration$preprocessing)) stopifnot(configuration$preprocessing == preprocessing)
if (preprocessing == "outer_fold") {
  stopifnot(identical(fit$tate_fit$preprocessing$mode, "outer_fold"),
    length(fit$tate_fit$preprocessing$parameters) == metadata$n_folds)
}
stopifnot(unname(metadata$target_site) == "Private", metadata$nuisance_lambda_rule == "min",
  metadata$crossfit_layers == 2L, metadata$aggregation_mode == "joint_tate",
  metadata$aggregation_cutoff == 2, is.null(metadata$calibration_control$source_lambda_rules))
source_radius <- metadata$M_tau
target_radius <- metadata$calibration_control$target_radius
if (is.null(target_radius)) target_radius <- source_radius
stopifnot(is.finite(source_radius), source_radius > 0, metadata$M_tau_inference == source_radius,
  is.finite(target_radius), target_radius > 0)
if (!is.null(configuration$source_radius)) stopifnot(configuration$source_radius == source_radius)
if (!is.null(configuration$target_radius)) stopifnot(configuration$target_radius == target_radius)
prefix <- configuration$asset_prefix
if (is.null(prefix)) prefix <- "rhc_private_min"
stopifnot(length(prefix) == 1L, grepl("^[a-z0-9_]+$", prefix))
asset <- function(suffix) file.path(output, paste0(prefix, suffix))
raw <- load_rhc_raw()
excluded <- c("No insurance", "Medicare & Medicaid")
cohort <- build_rhc_cohort(raw, outcome = "death30", site_var = "ninsclas",
  site_recode = setNames(rep(NA_character_, length(excluded)), excluded))
data <- build_rhc_data_split(cohort, K = 4L, target_site = "Private", seed = 42L)
stopifnot(identical(digest::digest(data, algo = "sha256"), configuration$data_sha256),
  nrow(cohort) == 5039L, all(vapply(data, function(site) ncol(site$W_outcome) == 61L &&
    identical(site$W_outcome, site$Z_site), logical(1L))))
methods <- read.csv(file.path(selection, "methods.csv"), stringsAsFactors = FALSE)
method_order <- c("RoCE", "Target-only", "Calibrated target-only", "sample_size", "inverse_variance", "federated_dr", "pooled_dr")
stopifnot(nrow(methods) == 7L, setequal(methods$method, method_order),
  all(is.finite(as.matrix(methods[-1L]))), all(methods$se > 0),
  max(abs(methods$ci_lower - methods$estimate + qnorm(.975) * methods$se)) < 1e-12,
  max(abs(methods$ci_upper - methods$estimate - qnorm(.975) * methods$se)) < 1e-12)
methods <- methods[match(method_order, methods$method), ]
stopifnot(abs(methods$estimate[1L] - fit$tate_fit$estimate) < 1e-12,
          abs(methods$se[1L] - fit$tate_fit$se) < 1e-12)
arms <- read.csv(file.path(selection, "arm_summary.csv"), stringsAsFactors = FALSE)
stopifnot(abs(arms$estimate[arms$arm == "mu1"] - arms$estimate[arms$arm == "mu0"] - methods$estimate[1L]) < 1e-12,
  abs(sum(arms$variance) - 2 * arms$covariance_mu1_mu0[1L] - methods$se[1L]^2) < 1e-12)
source_weights <- read.csv(file.path(selection, "sources.csv"), stringsAsFactors = FALSE)
mean_weights <- colMeans(fit$tate_fit$fold_weights)
for (arm in c("mu1", "mu0")) stopifnot(max(abs(source_weights[[paste0("weight_", arm)]] -
  mean_weights[paste(arm, source_weights$source, sep = ":")])) < 1e-12)
mapping <- attr(data, "site_mapping")
sites <- do.call(rbind, lapply(names(data), function(site) {
  sample <- data[[site]]
  index <- if (site == "t") "target" else site
  weights <- if (site == "t") c(1-sum(source_weights$weight_mu1), 1-sum(source_weights$weight_mu0)) else
    as.numeric(source_weights[source_weights$source == site, c("weight_mu1", "weight_mu0")])
  data.frame(site = site, insurance = unname(mapping[index]), n = sample$n,
    treated_n = sum(sample$A), deaths_n = sum(sample$Y), rhc_proportion = mean(sample$A),
    mortality_proportion = mean(sample$Y), weight_mu1 = weights[1L], weight_mu0 = weights[2L])
}))
stopifnot(sum(sites$n) == 5039L, max(abs(colSums(sites[c("weight_mu1", "weight_mu0")]) - 1)) < 1e-12)
dir.create(output, recursive = TRUE)
# Keep the calibrated anchor as the single displayed Target-only comparator.
# Original method identities remain explicit in the export provenance.
display_sources <- c("Target-only" = "Calibrated target-only", "SS" = "sample_size",
  "IVW" = "inverse_variance", "Federated-DR" = "federated_dr", "Pooled-DR" = "pooled_dr", "RoCE" = "RoCE")
methods <- methods[match(unname(display_sources), methods$method), ]
methods$method <- names(display_sources)
anchor <- fit$methods[as.character(fit$methods$method) == "Calibrated target-only", ]
stopifnot(nrow(methods) == 6L, nrow(anchor) == 1L,
  abs(methods$estimate[methods$method == "Target-only"] - anchor$estimate) < 1e-12,
  abs(methods$se[methods$method == "Target-only"] - anchor$se) < 1e-12)
write.csv(methods, asset("_tate.csv"), row.names = FALSE)
write.csv(sites, asset("_sites.csv"), row.names = FALSE)
write.csv(arms, asset("_arms.csv"), row.names = FALSE)
exclusions <- data.frame(insurance = excluded, n = vapply(excluded,
  function(label) sum(raw$ninsclas == label), integer(1L)))
write.csv(exclusions, file.path(output, "rhc_cohort_exclusions.csv"), row.names = FALSE)
jsonlite::write_json(list(configuration = configuration,
  inputs_sha256 = c(provenance, setNames(list(digest::digest(file.path(primary, "checkpoints/analysis.rds"),
    file = TRUE, algo = "sha256")), file.path(primary, "checkpoints/analysis.rds"))),
  displayed_method_sources = as.list(display_sources),
  checks = list(data_hash = TRUE, seven_saved_methods_verified = TRUE, six_displayed_methods = TRUE,
    target_only_is_calibrated_anchor = TRUE, normal_intervals = TRUE,
    arm_contrast = TRUE, cross_arm_variance = TRUE, mean_fold_weights = TRUE),
  source_radius = source_radius, target_radius = target_radius,
  preprocessing = preprocessing,
  scope = sprintf("Private target, min rule, source radius %g, target radius %g, %s preprocessing; fixed original partition and matched baselines",
    source_radius, target_radius, preprocessing)),
  asset("_provenance.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")

# Use the same method order, colours and typography as the simulation figures.
suppressPackageStartupMessages(library(ggplot2))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(script_path)), "publication_plot_style.R"))
stopifnot(identical(names(display_sources), publication_method_order))
values <- methods
values[-1L] <- 100 * values[-1L]
values$method <- factor(values$method, levels = rev(names(display_sources)))
figure <- ggplot(values, aes(y = method, x = estimate, colour = method, shape = method)) +
  geom_vline(xintercept = 0, colour = "#7D858D", linetype = "dashed", linewidth = .4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = .20, linewidth = .7) +
  geom_point(size = 2.7, stroke = .6) +
  scale_colour_manual(values = publication_colors, guide = "none") +
  scale_shape_manual(values = publication_shapes, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
  labs(x = "Mortality risk difference (percentage points)", y = NULL) +
  theme_publication()
ggsave(asset("_tate.pdf"), figure, width = 6.5, height = 3.2,
       device = grDevices::pdf, useDingbats = FALSE)
writeLines("PRIVATE_MIN_COHORT_METHODS_INTERVALS_AND_WEIGHTS_VERIFIED", file.path(output, "CHECKS_PASSED"))
print(sites)
print(exclusions)
print(methods)

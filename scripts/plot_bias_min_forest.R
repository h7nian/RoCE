# scripts/plot_bias_min_forest.R
#
# Regenerates the forest plot for the bias-minimum FACE-HD real-data
# configuration on the RHC cohort:
#   site_var = "ninsclas", K = 5, M_tau_inference = 5, extended covariates
#
# Reads the persisted result list from results/real_data/ and writes a
# paper-ready PDF with an explicit "bias-minimum configuration" subtitle
# and JASA-friendly sizing (7" x 5").

suppressPackageStartupMessages({
  if (!requireNamespace("FACEHD", quietly = TRUE)) {
    stop("FACEHD package is not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for forest plot regeneration.")
  }
  library(FACEHD)
})

setting_id   <- "rhc_K5_death30_kf10_A1_ninsclas_Mt5"
results_dir  <- "results/real_data"
rds_path     <- file.path(results_dir, paste0(setting_id, ".rds"))
pdf_out      <- file.path(results_dir, paste0(setting_id, "_forest_paper.pdf"))

if (!file.exists(rds_path)) {
  stop("Missing result RDS: ", rds_path)
}

result <- readRDS(rds_path)

target_only_val <- result$methods$estimate[result$methods$method == "target_only"]
facehd_val      <- result$methods$estimate[result$methods$method == "facehd"]
bias_val        <- facehd_val - target_only_val

forest <- plot_forest_methods(
  result$methods,
  title    = "RHC real-data application: 30-day mortality under RHC exposure",
  subtitle = sprintf(
    "Bias-minimum FACE-HD configuration (site_var = ninsclas, K = 5, M_tau = 5, p ~ 75); bias vs. target-only = %+.3f",
    bias_val
  ),
  xlab     = "Estimate of target-site potential-outcome mean mu^1_t"
)

save_plot(forest, pdf_out, width = 7, height = 5)

cat(sprintf("[OK] Forest plot written to %s\n", pdf_out))
cat(sprintf("    target_only = %.4f, facehd = %.4f, bias = %+.4f\n",
            target_only_val, facehd_val, bias_val))

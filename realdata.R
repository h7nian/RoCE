# realdata.R - RoCE Real-Data Experiment (Right-Heart Catheterization)
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# Potential-outcome mean at the target site:
#   mu^A_t = E_t[Y(A)] for A in {0, 1}
# with A = 1 corresponding to right-heart catheterization (RHC exposure).
#
# Pipeline (parallel in spirit to main.R):
#   1. Parse args.
#   2. Load installed RoCE package (devtools::load_all fallback for dev).
#   3. Run the end-to-end RHC experiment via run_rhc_experiment().
#   4. Persist artefacts under results/real_data/:
#        - <setting>.rds              full result list
#        - <setting>_methods.csv      method-vs-estimate table
#        - <setting>_forest.pdf       forest plot
#        - <setting>_weights.pdf      aggregation-weight bar chart
#        - <setting>_pairwise.pdf     pairwise-vs-target scatter
#      Plots are skipped gracefully when ggplot2 is not available.
# =============================================================================

args <- commandArgs(TRUE)
print("Command line arguments:")
print(args)

# ============================================================================
# Command-Line Arguments (positional; all optional after the first)
# ============================================================================
# Order puts always-provided arguments first and keeps the optional empty-
# eligible arguments at the very end so that the shell's word-splitting of
# "${TARGET_SITE:-}" (possibly empty) cannot swallow a later value.
#
# arg1: K                number of sites (default 5)
# arg2: outcome          "death30" (default), "death180", or "los"
# arg3: n_folds          outer cross-fitting folds (default 10)
# arg4: A_val            target treatment level (0 or 1, default 1)
# arg5: job_id           SLURM job id (always passed; falls back to timestamp)
# arg6: site_var         raw RHC column for site partitioning (default "cat1")
# arg7: M_tau_inference  inference-time truncation cap on phi^T gamma
#                        (default 3; accepts "Inf" for no truncation)
# arg8: target_site      site-level value overriding the default target
#                        (optional; MUST be last so empty expansion cannot
#                        shift earlier indices)
# ============================================================================
.parse_m_tau <- function(x) {
  if (is.null(x) || !nzchar(x)) return(3.0)
  if (tolower(x) %in% c("inf", "infinity")) return(Inf)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v) || v <= 0) stop(sprintf("Invalid M_tau_inference = '%s'", x))
  v
}

K_arg             <- if (length(args) >= 1) as.integer(args[1]) else 5L
outcome_arg       <- if (length(args) >= 2) tolower(args[2]) else "death30"
n_folds_arg       <- if (length(args) >= 3) as.integer(args[3]) else 10L
A_val_arg         <- if (length(args) >= 4) as.integer(args[4]) else 1L
job_id            <- if (length(args) >= 5 && nzchar(args[5])) args[5] else format(Sys.time(), "%Y%m%d_%H%M%S")
site_var_arg      <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "cat1"
m_tau_arg         <- .parse_m_tau(if (length(args) >= 7) args[7] else NULL)
target_site_arg   <- if (length(args) >= 8 && nzchar(args[8])) args[8] else NULL

stopifnot(
  "K must be an integer >= 2"            = !is.na(K_arg) && K_arg >= 2L,
  "n_folds must be an integer >= 3"      = !is.na(n_folds_arg) && n_folds_arg >= 3L,
  "A_val must be 0 or 1"                 = A_val_arg %in% c(0L, 1L),
  "outcome must be one of death30/death180/los" = outcome_arg %in% c("death30", "death180", "los"),
  "site_var must be a non-empty character" = is.character(site_var_arg) && nzchar(site_var_arg)
)

# Translate the SITE_PRESET environment variable into a named character
# vector consumed by build_rhc_cohort(site_recode = ...). Presets are
# intentionally hardcoded here (rather than introducing a dynamic parser)
# so that each named preset is explicit, documented, and reviewable in
# source control.
site_preset_arg <- Sys.getenv("SITE_PRESET", "")
site_recode_arg <- if (!nzchar(site_preset_arg)) {
  NULL
} else {
  switch(
    site_preset_arg,
    "ninsclas4" = c(
      "No insurance"        = NA_character_,
      "Medicaid"            = "Medicaid_any",
      "Medicare & Medicaid" = "Medicaid_any"
    ),
    stop(sprintf("Unknown SITE_PRESET = '%s' (known: ninsclas4)", site_preset_arg))
  )
}

# Map the PHI environment variable to an actual basis-expansion function.
# Hardcoded presets (rather than user-supplied code) keep the command-line
# surface declarative and reviewable.
#
# RoCE is loaded below (via library(RoCE)), so we cannot yet reference
# phi_rhc_bspline by its bare name. Use the fully-qualified namespace
# accessor instead; this also makes the dependency on RoCE explicit.
phi_preset_arg <- Sys.getenv("PHI", "identity")
phi_fun <- switch(
  phi_preset_arg,
  "identity" = base::identity,
  "bspline"  = RoCE::phi_rhc_bspline,
  stop(sprintf("Unknown PHI = '%s' (known: identity, bspline)", phi_preset_arg))
)

# ============================================================================
# Load RoCE Package (mirrors main.R)
# ============================================================================
# quietly = FALSE so the actual loadNamespace error (e.g. ABI mismatch from a
# wrong R module version) surfaces in the log instead of silently falling to
# the misleading "neither found" branch.
if (requireNamespace("RoCE", quietly = FALSE)) {
  library(RoCE)
} else if (requireNamespace("devtools", quietly = FALSE)) {
  cat("RoCE not installed; loading via devtools::load_all()...\n")
  devtools::load_all(".")
} else {
  stop("Neither installed RoCE package nor devtools found. ",
       "Install the package with: R CMD INSTALL .")
}

# ============================================================================
# Setting Identifier and Output Paths
# ============================================================================
# M_tau suffix: "3" -> "Mt3", "10" -> "Mt10", Inf -> "Mtinf".
.m_tau_tag <- function(m) {
  if (is.infinite(m)) return("Mtinf")
  # Drop trailing ".0" for whole numbers; otherwise keep one decimal place.
  if (m == as.integer(m)) sprintf("Mt%d", as.integer(m)) else sprintf("Mt%.1f", m)
}

setting_id <- sprintf("rhc_K%d_%s_kf%d_A%d_%s_%s",
                      K_arg, outcome_arg, n_folds_arg, A_val_arg,
                      site_var_arg, .m_tau_tag(m_tau_arg))
if (nzchar(site_preset_arg)) {
  setting_id <- paste0(setting_id, "_", site_preset_arg)
}
if (phi_preset_arg != "identity") {
  setting_id <- paste0(setting_id, "_", phi_preset_arg)
}

results_dir <- file.path("results", "real_data")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

rds_path      <- file.path(results_dir, paste0(setting_id, ".rds"))
methods_csv   <- file.path(results_dir, paste0(setting_id, "_methods.csv"))
pairwise_csv  <- file.path(results_dir, paste0(setting_id, "_pairwise.csv"))
weights_csv   <- file.path(results_dir, paste0(setting_id, "_weights.csv"))
forest_pdf    <- file.path(results_dir, paste0(setting_id, "_forest.pdf"))
weights_pdf   <- file.path(results_dir, paste0(setting_id, "_weights.pdf"))
pairwise_pdf  <- file.path(results_dir, paste0(setting_id, "_pairwise.pdf"))

cat(sprintf("==============================================\n"))
cat(sprintf("RoCE Real-Data Experiment: RHC\n"))
cat(sprintf("==============================================\n"))
cat(sprintf("Job ID:        %s\n", job_id))
cat(sprintf("Setting ID:    %s\n", setting_id))
cat(sprintf("K:             %d\n", K_arg))
cat(sprintf("Outcome:       %s\n", outcome_arg))
cat(sprintf("n_folds:       %d\n", n_folds_arg))
cat(sprintf("A_val:         %d\n", A_val_arg))
cat(sprintf("Site var:      %s\n", site_var_arg))
cat(sprintf("M_tau_inf:     %s\n", format(m_tau_arg)))
cat(sprintf("Site preset:   %s\n",
            if (nzchar(site_preset_arg)) site_preset_arg else "<none>"))
cat(sprintf("Phi basis:     %s\n", phi_preset_arg))
cat(sprintf("Target site:   %s\n",
            if (is.null(target_site_arg))
              sprintf("<auto, largest level of %s>", site_var_arg)
            else target_site_arg))
cat(sprintf("==============================================\n\n"))

# ============================================================================
# Run Experiment
# ============================================================================
result <- run_rhc_experiment(
  K               = K_arg,
  target_site     = target_site_arg,
  outcome         = outcome_arg,
  site_var        = site_var_arg,
  site_recode     = site_recode_arg,
  n_folds         = n_folds_arg,
  A_val           = A_val_arg,
  phi             = phi_fun,
  n_cores         = -1L,        # use all-but-one core on the SLURM node
  M_tau_inference = m_tau_arg,
  verbose         = TRUE
)

# ============================================================================
# Persist Tabular Artefacts
# ============================================================================
saveRDS(result, file = rds_path)
cat(sprintf("[OK] Result list written to %s\n", rds_path))

utils::write.csv(result$methods, file = methods_csv, row.names = FALSE)
utils::write.csv(result$pairwise, file = pairwise_csv, row.names = FALSE)

# Source labels come from the pairwise table; this is the canonical label
# vector used throughout the rest of the pipeline. It matches names(result$weights)
# when weights are named, and falls back cleanly when they are not.
source_labels <- result$pairwise$source
weights_df <- data.frame(
  source        = source_labels,
  weight        = as.numeric(result$weights),
  pairwise_d_sq = result$metadata$pairwise_d_sq,
  stringsAsFactors = FALSE
)
utils::write.csv(weights_df, file = weights_csv, row.names = FALSE)
cat(sprintf("[OK] CSV tables written to %s\n", results_dir))

# Summary printed to log so SLURM stdout captures the numeric headline.
cat("\n=== Methods table ===\n")
print(result$methods, digits = 4, row.names = FALSE)
cat("\n=== RoCE aggregation weights ===\n")
print(weights_df, digits = 4, row.names = FALSE)

# ============================================================================
# Plots (optional; skipped without ggplot2)
# ============================================================================
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
if (!has_ggplot2) {
  cat("\n[INFO] ggplot2 not available; skipping PDF figures.\n",
      "       Tabular outputs above contain all information needed for later plotting.\n",
      sep = "")
} else {
  forest <- plot_forest_methods(
    result$methods,
    title    = NULL,
    subtitle = NULL,
    xlab     = "Potential-outcome mean"
  )
  save_plot(forest, forest_pdf, width = 7, height = 4.5)
  cat(sprintf("[OK] Forest plot: %s\n", forest_pdf))

  weights_plot <- plot_aggregation_weights(
    weights      = result$weights,
    source_labels = source_labels,
    penalty_d2   = result$metadata$pairwise_d_sq,
    title        = "RoCE aggregation weights",
    subtitle     = sprintf("Zero weights = selection penalty dominated by (%s gap)^2", outcome_arg)
  )
  save_plot(weights_plot, weights_pdf, width = 7, height = 4.0)
  cat(sprintf("[OK] Weights plot: %s\n", weights_pdf))

  pairwise_plot <- plot_pairwise_vs_target(
    target_only_estimate = result$target_only,
    pairwise             = transform(result$pairwise, source = source_labels),
    title                = "Pairwise source-assisted vs. target-only",
    subtitle             = "Diagonal: source unbiased for target parameter"
  )
  save_plot(pairwise_plot, pairwise_pdf, width = 6, height = 5)
  cat(sprintf("[OK] Pairwise plot: %s\n", pairwise_pdf))
}

cat(sprintf("\n==============================================\n"))
cat(sprintf("Real-data job completed at: %s\n", format(Sys.time())))
cat(sprintf("==============================================\n"))

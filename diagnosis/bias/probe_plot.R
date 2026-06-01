#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/probe_plot.R
# ----------------------------------------------------------------------------
# Reads probe_scan_long_*.csv, probe_scan_spread_*.csv (and probe_kf_*.csv
# if present) and produces base-R PNG plots that visualise the fold-to-fold
# nuisance instability we have been quantifying.
#
# Plots:
#   1. mu_pred_ts vs k2, faceted by (n, p) — fold-to-fold instability
#   2. mu_pred_ts_spread vs p / n_calib_arm (log x) — the "cliff" scatter
#   3. Explosion rate (spread > 0.2) bar chart per setting
#   4. (If probe_kf data) spread vs K_f for catastrophic/production settings
#
# Output: diagnosis/bias/probe_out/plots/*.png
# ============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

out_dir  <- "diagnosis/bias/probe_out"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Find newest probe_scan_* and probe_kf_* files
newest <- function(pattern) {
  fs <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  if (length(fs) == 0L) return(NULL)
  fs[which.max(file.mtime(fs))]
}
scan_long_f   <- newest("^probe_scan_long_.*\\.csv$")
scan_spread_f <- newest("^probe_scan_spread_.*\\.csv$")
kf_long_f     <- newest("^probe_kf_long_.*\\.csv$")
kf_spread_f   <- newest("^probe_kf_spread_.*\\.csv$")

cat(sprintf("[plot] scan_long   = %s\n",  scan_long_f   %||% "(missing)"))
cat(sprintf("[plot] scan_spread = %s\n",  scan_spread_f %||% "(missing)"))
cat(sprintf("[plot] kf_long     = %s\n",  kf_long_f     %||% "(missing)"))
cat(sprintf("[plot] kf_spread   = %s\n",  kf_spread_f   %||% "(missing)"))

# ============================================================================
# Plot 1: mu_pred_ts vs k2, faceted by (n, p)
# ============================================================================
if (!is.null(scan_long_f)) {
  long <- read.csv(scan_long_f, stringsAsFactors = FALSE)
  long$setting <- sprintf("n=%d, p=%d", long$n, long$p)
  settings <- unique(long[order(long$n, long$p), c("setting", "n", "p")])

  pdf(file.path(plot_dir, "1_mu_pred_ts_vs_k2.pdf"),
      width = 12, height = 8)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (ii in seq_len(nrow(settings))) {
    sub <- long[long$setting == settings$setting[ii], ]
    plot(NA, NA, xlim = range(sub$k2), ylim = c(0, 1),
         xlab = "k2 (secondary fold)", ylab = "mu_pred_ts (single source × k1=1)",
         main = settings$setting[ii], cex.main = 1.1)
    abline(h = 0.5, lty = 3, col = "gray60")
    # one line per (seed, source)
    grp <- interaction(sub$seed, sub$source, drop = TRUE)
    for (g in unique(grp)) {
      ss <- sub[grp == g, ]
      ss <- ss[order(ss$k2), ]
      lines(ss$k2, ss$mu_pred_ts,
            col = adjustcolor("steelblue", alpha.f = 0.5), lwd = 1.5)
      points(ss$k2, ss$mu_pred_ts,
             pch = 16, col = adjustcolor("steelblue", alpha.f = 0.5))
    }
    legend("topright", legend = sprintf("3 seeds × %d sources", length(unique(sub$source))),
           bty = "n", cex = 0.85)
  }
  mtext("mu_pred_ts across secondary folds (k2 sweep, k1=1)", outer = TRUE,
        cex = 1.2, font = 2)
  dev.off()
  cat("[plot] wrote 1_mu_pred_ts_vs_k2.pdf\n")
}

# ============================================================================
# Plot 2: spread vs p/n_calib_arm (the "cliff")
# ============================================================================
if (!is.null(scan_spread_f)) {
  spread <- read.csv(scan_spread_f, stringsAsFactors = FALSE)
  spread$setting <- sprintf("n=%d, p=%d", spread$n, spread$p)

  pdf(file.path(plot_dir, "2_spread_vs_ratio.pdf"),
      width = 8, height = 6)
  par(mar = c(5, 5, 3, 1))
  cols <- rainbow(length(unique(spread$setting)))
  ucol <- setNames(cols, unique(spread$setting))
  plot(spread$ratio_p_over_calib_arm, spread$mu_pred_ts_spread,
       log = "x",
       xlab = "p / n_calib_arm  (log scale)",
       ylab = "mu_pred_ts spread (max - min across k2)",
       main = "Nuisance fold-to-fold spread vs sample-budget ratio",
       pch = 16, cex = 1.4,
       col = adjustcolor(ucol[spread$setting], alpha.f = 0.7))
  abline(h = 0.2, lty = 2, col = "red")
  abline(v = 0.2, lty = 2, col = "darkgreen")
  legend("topleft", legend = names(ucol), col = ucol, pch = 16,
         bty = "n", cex = 0.9)
  text(0.21, 0.95, "ratio = 0.2 (rough safe boundary)", srt = 90,
       col = "darkgreen", adj = c(0, -0.5), cex = 0.8)
  text(grconvertX(0.95, "npc"), 0.22, "spread = 0.2 (blow-up threshold)",
       col = "red", adj = c(1, -0.5), cex = 0.8)
  dev.off()
  cat("[plot] wrote 2_spread_vs_ratio.pdf\n")
}

# ============================================================================
# Plot 3: explosion rate bar chart
# ============================================================================
if (!is.null(scan_spread_f)) {
  spread <- read.csv(scan_spread_f, stringsAsFactors = FALSE)
  spread$is_blown <- spread$mu_pred_ts_spread > 0.2
  rate <- aggregate(is_blown ~ tag + n + p, data = spread, FUN = mean)
  rate <- rate[order(rate$n, rate$p), ]
  rate$label <- sprintf("%s\nn=%d,p=%d", rate$tag, rate$n, rate$p)

  pdf(file.path(plot_dir, "3_explosion_rate.pdf"),
      width = 8, height = 5)
  par(mar = c(7, 5, 3, 1))
  bp <- barplot(rate$is_blown * 100,
                names.arg = rate$label, las = 2,
                ylim = c(0, 110),
                ylab = "% (seed, source) pairs with spread > 0.2",
                main = "Explosion rate by (n, p) setting",
                col = ifelse(rate$is_blown == 1, "firebrick",
                             ifelse(rate$is_blown > 0.5, "orange", "steelblue")),
                border = NA, cex.names = 0.85)
  text(bp, rate$is_blown * 100 + 3,
       labels = sprintf("%.0f%%", rate$is_blown * 100), cex = 0.9)
  dev.off()
  cat("[plot] wrote 3_explosion_rate.pdf\n")
}

# ============================================================================
# Plot 4: K_f sweep (only if probe_kf data is present)
# ============================================================================
if (!is.null(kf_spread_f)) {
  kfs <- read.csv(kf_spread_f, stringsAsFactors = FALSE)

  pdf(file.path(plot_dir, "4_kf_sweep.pdf"),
      width = 10, height = 5)
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), oma = c(0, 0, 2, 0))
  for (st in sort(unique(kfs$setting_tag))) {
    sub <- kfs[kfs$setting_tag == st, ]
    agg <- aggregate(mu_pred_ts_spread ~ K_f, data = sub, FUN = median)
    rng <- range(sub$mu_pred_ts_spread, na.rm = TRUE)
    plot(agg$K_f, agg$mu_pred_ts_spread,
         type = "b", pch = 16, cex = 1.5, lwd = 2, col = "steelblue",
         xlab = "K_f (number of cross-fitting folds)",
         ylab = "median mu_pred_ts spread",
         main = st, ylim = c(0, max(rng, 1) * 1.05),
         xaxt = "n")
    axis(1, at = sort(unique(sub$K_f)))
    abline(h = 0.2, lty = 2, col = "red")
    # also overlay individual (seed, source) points
    points(sub$K_f, sub$mu_pred_ts_spread,
           pch = 21, col = adjustcolor("gray40", alpha.f = 0.5),
           bg = adjustcolor("white", alpha.f = 0.3))
  }
  mtext("Effect of K_f on fold-to-fold spread", outer = TRUE,
        cex = 1.2, font = 2)
  dev.off()
  cat("[plot] wrote 4_kf_sweep.pdf\n")
} else {
  cat("[plot] kf_spread file missing — skipping plot 4 (re-run after probe_kf finishes)\n")
}

cat("\n[plot] all plots in: ", plot_dir, "\n", sep = "")

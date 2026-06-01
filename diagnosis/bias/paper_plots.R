#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/paper_plots.R
# ----------------------------------------------------------------------------
# Publication-quality figures (PDF, vector) from the FACE-HD simulation
# summaries. Pulls from BOTH:
#   - results/*_summary.csv        (production grid; K_f parsed from filename)
#   - diagnosis/bias/verify_out/*  (K_f sweep; n_folds column present)
#
# Figures (results/figures/*.pdf) — K_f=3 results only (no K_f=10 baseline):
#   fig1_coverage_by_method   coverage of all methods × config at the
#                             production setting (n=5000, K=3, p=10, K_f=3)
#   fig2_coverage_vs_dim      one/two-round coverage vs p (dimension), K_f=3
#
# MSI's R has no png/cairo (capabilities() all FALSE) so output is PDF —
# which is the preferred vector format for papers anyway.
# ============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

fig_dir <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Method display labels + ordering + styling ---------------------------
METHOD_LABEL <- c(
  one_round_crossfit = "FACE-HD (1-round)",
  two_round_crossfit = "FACE-HD (2-round)",
  federated_dr       = "Fed-DR",
  pooled_dr          = "Pooled-DR",
  oracle_dr          = "Oracle-DR",
  target_only        = "Target-only",
  tilted_aipw        = "Tilted-AIPW",
  inverse_variance   = "IVW",
  sample_size        = "Size-weighted"
)
# our two methods get emphasis
OUR_METHODS <- c("one_round_crossfit", "two_round_crossfit")
METHOD_ORDER <- c("one_round_crossfit", "two_round_crossfit",
                  "federated_dr", "pooled_dr", "oracle_dr",
                  "tilted_aipw", "inverse_variance", "sample_size",
                  "target_only")

# ---- Load + tag production summaries --------------------------------------
parse_production <- function() {
  fs <- list.files("results", pattern = "_summary\\.csv$", full.names = TRUE)
  fs <- fs[grepl("^n[0-9]+_K[0-9]+_p[0-9]+_C[0-9]", basename(fs))]
  if (length(fs) == 0) return(NULL)
  rows <- lapply(fs, function(f) {
    b <- basename(f)
    m <- regmatches(b, regexec(
      "^n([0-9]+)_K([0-9]+)_p([0-9]+)_(C[0-9])_.*_kf([0-9]+)_", b))[[1]]
    if (length(m) < 6) return(NULL)
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$n_total <- as.integer(m[2]); df$K <- as.integer(m[3])
    df$p <- as.integer(m[4]); df$config <- m[5]
    df$n_folds <- as.integer(m[6])
    df$source_file <- b
    df
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ---- Load + tag verify_out summaries (K_f sweep) --------------------------
parse_verify <- function() {
  fs <- list.files("diagnosis/bias/verify_out", pattern = "_summary\\.csv$",
                   full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  rows <- lapply(fs, function(f) {
    b <- basename(f)
    m <- regmatches(b, regexec(
      "verify_(C[0-9])_n([0-9]+)_K([0-9]+)_p([0-9]+)_kf([0-9]+)_", b))[[1]]
    if (length(m) < 6) return(NULL)
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$config <- m[2]; df$n_total <- as.integer(m[3]); df$K <- as.integer(m[4])
    df$p <- as.integer(m[5]); df$n_folds <- as.integer(m[6])
    df$source_file <- b
    df[, intersect(names(df), c("method","config","n_total","K","p","n_folds",
                                "bias_mean","coverage","source_file"))]
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

prod <- parse_production()
ver  <- parse_verify()
cat(sprintf("[paper_plots] production rows: %d ; verify rows: %d\n",
            nrow(prod) %||% 0, nrow(ver) %||% 0))

mlab <- function(m) METHOD_LABEL[m] %||% m
is_ours <- function(m) m %in% OUR_METHODS

# ============================================================================
# FIG 1 — coverage by method × config (n=5000, K=3, K_f=3), one figure per p
# ============================================================================
if (!is.null(prod)) {
  base_all <- prod[prod$n_total == 5000 & prod$K == 3 & prod$n_folds == 3 &
                   prod$config %in% c("C1","C2","C3"), ]
  p_values <- sort(unique(base_all$p))
  for (pp in p_values) {
    base <- base_all[base_all$p == pp, ]
    if (nrow(base) == 0) next
    configs <- sort(unique(base$config))
    out_pdf <- file.path(fig_dir, sprintf("fig1_coverage_by_method_p%d.pdf", pp))
    pdf(out_pdf, width = 11, height = 4.2)
    par(mfrow = c(1, length(configs)), mar = c(8, 4.2, 3, 1), oma = c(0,0,2,0))
    for (cfg in configs) {
      sub <- base[base$config == cfg, ]
      meths <- METHOD_ORDER[METHOD_ORDER %in% sub$method]
      kf3  <- sapply(meths, function(m){ v<-sub$coverage[sub$method==m]; if(length(v)) v[1] else NA})
      x <- seq_along(meths)
      plot(NA, NA, xlim = c(0.5, length(meths)+0.5), ylim = c(0.80, 1.0),
           xaxt = "n", xlab = "", ylab = "Coverage", main = cfg, cex.main = 1.2)
      abline(h = 0.95, lty = 2, col = "gray40")
      rect(0.5, 0.94, length(meths)+0.5, 0.96, col = adjustcolor("gray80", 0.3), border = NA)
      cols <- ifelse(is_ours(meths), "firebrick", "steelblue4")
      points(x, kf3, pch = 19, cex = 1.7, col = cols)
      axis(1, at = x, labels = mlab(meths), las = 2, cex.axis = 0.8)
    }
    mtext(sprintf("Coverage by method (n=5000, K=3, p=%d, K_f=3).  Red = FACE-HD.  Band = [0.94, 0.96].", pp),
          outer = TRUE, cex = 0.95, font = 2)
    dev.off()
    cat(sprintf("[paper_plots] wrote %s\n", basename(out_pdf)))
  }
}

# ============================================================================
# FIG 2 (coverage vs dimension p, K_f=3, our methods), faceted by config
# (renumbered: K_f-ablation and bias-shrink figures dropped — they were
#  K_f=10-vs-K_f=3 comparisons and the K_f=10 baseline is not plotted.)
# ============================================================================
if (!is.null(prod)) {
  d3 <- prod[prod$n_folds == 3 & prod$n_total == 5000 & prod$K == 3 &
             prod$method %in% OUR_METHODS, ]
  if (nrow(d3) > 0 && length(unique(d3$p)) >= 2) {
    configs <- sort(unique(d3$config))
    pdf(file.path(fig_dir, "fig2_coverage_vs_dim.pdf"),
        width = 4 * length(configs), height = 4.2)
    par(mfrow = c(1, length(configs)), mar = c(4.5, 4.5, 3, 1))
    cols <- c(one_round_crossfit = "firebrick", two_round_crossfit = "navy")
    pchs <- c(one_round_crossfit = 19, two_round_crossfit = 17)
    ps <- sort(unique(d3$p))
    for (cfg in configs) {
      plot(NA, NA, xlim = range(ps), ylim = c(0.85, 1.0), log = "x",
           xlab = "p (covariate dimension)", ylab = "Coverage",
           main = sprintf("%s (K_f=3)", cfg), xaxt = "n")
      axis(1, at = ps)
      abline(h = 0.95, lty = 2, col = "gray40")
      for (m in OUR_METHODS) {
        sm <- d3[d3$config == cfg & d3$method == m, ]
        sm <- sm[order(sm$p), ]
        if (nrow(sm)) {
          lines(sm$p, sm$coverage, col = cols[m], lwd = 2)
          points(sm$p, sm$coverage, col = cols[m], pch = pchs[m], cex = 1.6)
        }
      }
      if (cfg == configs[1])
        legend("bottomleft", legend = mlab(OUR_METHODS),
               col = cols[OUR_METHODS], pch = pchs[OUR_METHODS], lwd = 2, bty = "n", cex = 0.85)
    }
    dev.off()
    cat("[paper_plots] wrote fig2_coverage_vs_dim.pdf\n")
  } else {
    cat("[paper_plots] fig2 (vs dim) skipped (need >=2 p values at K_f=3; have:",
        paste(sort(unique(d3$p)), collapse=","), ")\n")
  }
}

cat("\n[paper_plots] figures in: ", fig_dir, "\n", sep = "")
cat("[paper_plots] DONE\n")

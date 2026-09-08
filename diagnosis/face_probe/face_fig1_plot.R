#!/usr/bin/env Rscript
# FACE.pdf Figure-1 analog: Bias / RMSE / Coverage / 95%-CI length versus the
# source-deviation gradient rho, faceted by the number of source sites K.
# `direct` mode requires TATE rows for every method. `legacy` mode requires the
# archived treated-potential-outcome-mean rows for every method. Keeping the
# method maps separate prevents accidental mixed-estimand plots.
suppressPackageStartupMessages(library(ggplot2))

val_dir <- Sys.getenv(
  "ROCE_RESULTS_DIR",
  "results/direct_tate_mc500_b5000/summary"
)
fig_p <- Sys.getenv("ROCE_FIG_P", "")  # "" = all p mixed; e.g. "50" restricts to p=50
estimand_scope <- match.arg(
  Sys.getenv("ROCE_FIG_ESTIMAND_SCOPE", "direct"),
  c("direct", "legacy")
)
pat <- if (nzchar(fig_p)) sprintf("^negtransfer_.*_p%s_K[0-9]+\\.csv$", fig_p) else "^negtransfer_.*_K[0-9]+\\.csv$"
files <- list.files(val_dir, pattern = pat, full.names = TRUE)
if (length(files) == 0) stop("no negtransfer_*_K*.csv files found yet")
dat <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))

# Main Fig-1 is the well-specified config across source-site counts K; C2/C3 are
# kept for a separate misspecification panel (avoid conflating configs in one K facet).
main_cfg <- Sys.getenv("ROCE_FIG_CONFIG", "C1")
if ("config" %in% names(dat)) dat <- dat[dat$config == main_cfg, ]
if (nzchar(fig_p) && "p" %in% names(dat)) dat <- dat[dat$p == as.integer(fig_p), ]

meth_lv <- if (estimand_scope == "direct") {
  c("target_only_ate", "sample_size_ate", "inverse_variance_ate",
    "tilted_aipw_ate", "federated_dr_ate", "pooled_dr_ate",
    "one_round_crossfit_ate")
} else {
  # Archived arm-specific summaries use unsuffixed labels for every estimator,
  # including RoCE and target-only. Do not mix these with `_ate` TATE rows.
  c("target_only", "sample_size", "inverse_variance",
    "tilted_aipw", "federated_dr", "pooled_dr",
    "one_round_crossfit")
}
meth_lab <- c("Target-only","SS","IVW","Tilted-AIPW",
              "Federated-DR","Pooled-DR","RoCE")

# Tilted-AIPW fits all nuisances by UNPENALIZED MLE (lambda=0 density ratio); in
# high dimension this overfits, producing extreme tilting weights and a degenerate
# variance (mean CI ~1e4, coverage vacuously 1.0) inconsistent with the penalized
# DR baselines (federated_dr/pooled_dr use a CV lambda). It is informative only at
# low dimension, so we drop it from the high-dimensional (p>=50) figures.
if (nzchar(fig_p) && as.integer(fig_p) >= 50) {
  keep <- !grepl("^tilted_aipw(?:_ate)?$", meth_lv)
  meth_lv <- meth_lv[keep]; meth_lab <- meth_lab[keep]
}
missing_methods <- setdiff(meth_lv, unique(dat$method))
if (length(missing_methods) > 0L) {
  stop(
    estimand_scope, " summary is missing method(s): ",
    paste(missing_methods, collapse = ", "),
    if (estimand_scope == "direct") {
      ". Refusing to mix arm-level and TATE rows."
    } else {
      ". Refusing to relabel or alter archived legacy rows."
    }
  )
}
dat <- dat[dat$method %in% meth_lv, ]
dat$method <- factor(dat$method, levels = meth_lv, labels = meth_lab)
dat$Kf <- factor(sprintf("K = %d", dat$K), levels = sprintf("K = %d", sort(unique(dat$K))))

# Long format over the four FACE panels (base R; no tidyr dependency).
panels <- list(Bias = "bias_mean", RMSE = "rmse",
               Coverage = "coverage", `CI length` = "ci_width_mean")
long <- do.call(rbind, lapply(names(panels), function(nm) {
  data.frame(method = dat$method, rho = dat$rho, Kf = dat$Kf,
             metric = factor(nm, levels = names(panels)),
             value = dat[[panels[[nm]]]], stringsAsFactors = FALSE)
}))

# Each method gets a distinct colour AND a distinct line type (so the figure reads
# in black-and-white too); lines only, no point markers. Keyed by label so the
# mapping is stable whether or not Tilted-AIPW is present (low- vs high-dim).
col_map <- c("RoCE"="#D95F02", "Target-only"="#1B9E77", "SS"="#7570B3",
             "IVW"="#E7298A", "Federated-DR"="#66A61E", "Pooled-DR"="#E6AB02",
             "Tilted-AIPW"="#A6761D")
lt_map  <- c("RoCE"="solid", "Target-only"="longdash", "SS"="62",
             "IVW"="dashed", "Federated-DR"="dotdash", "Pooled-DR"="twodash",
             "Tilted-AIPW"="1232")

p <- ggplot(long, aes(rho, value, colour = method, linetype = method, group = method)) +
  geom_line(linewidth = 1.0) +
  facet_grid(metric ~ Kf, scales = "free_y") +
  scale_colour_manual(values = col_map, name = NULL) +
  scale_linetype_manual(values = lt_map, name = NULL) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  labs(x = expression("source deviation " * rho * " (log-odds)"), y = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 17),
        strip.text = element_text(size = 19),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        # Leave enough device padding for the large, rotated row-strip labels.
        # This avoids clipping the final label without sacrificing readability.
        plot.margin = margin(t = 12, r = 28, b = 20, l = 12, unit = "pt"))

out_pdf <- file.path(val_dir, if (nzchar(fig_p))
  sprintf("fig1_negtransfer_%s_p%s.pdf", main_cfg, fig_p) else "fig1_negtransfer.pdf")
ggsave(out_pdf, p, width = 3.2 + 2.5 * length(unique(dat$K)), height = 9.0)
cat(sprintf("wrote %s  (K = {%s}, %d methods, rho = {%s})\n", out_pdf,
            paste(sort(unique(dat$K)), collapse=","), nlevels(dat$method),
            paste(sort(unique(dat$rho)), collapse=",")))

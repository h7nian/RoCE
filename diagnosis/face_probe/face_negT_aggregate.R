#!/usr/bin/env Rscript
# Combine all chunk CSVs -> summarize per (config, K, rho) -> write the
# negtransfer_<cfg>_p<P>_K<K>.csv files consumed by face_fig1_plot.R.
suppressPackageStartupMessages(library(RoCE))
cdir <- "diagnosis/face_probe/validation/chunks"
files <- list.files(cdir, pattern = "^negT_K.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no chunk CSVs found")
raw <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
if (!"heterogeneity_type" %in% names(raw)) raw$heterogeneity_type <- "none"
cat(sprintf("aggregating %d chunk files, %d raw rows\n", length(files), nrow(raw)))

cells <- split(raw, list(raw$config, raw$K, raw$rho, raw$p), drop = TRUE)
summ <- do.call(rbind, lapply(cells, function(df) {
  s <- tryCatch(summarize_results(df), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  s$rho <- df$rho[1]; s$K <- df$K[1]; s$config <- df$config[1]; s$p <- df$p[1]
  s$n_chunksims <- length(unique(df$sim_id))
  s
}))
for (pp in sort(unique(summ$p)))
  for (cfg in unique(summ$config[summ$p == pp]))
    for (k in sort(unique(summ$K[summ$config == cfg & summ$p == pp]))) {
      sub <- summ[summ$config == cfg & summ$K == k & summ$p == pp, ]
      if (nrow(sub) == 0) next
      out <- sprintf("diagnosis/face_probe/validation/negtransfer_%s_p%d_K%d.csv", cfg, pp, k)
      write.csv(sub, out, row.names = FALSE)
      cat(sprintf("wrote %s (%d rows, sims/cell~%d)\n", out, nrow(sub), max(sub$n_chunksims)))
    }

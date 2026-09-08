#!/usr/bin/env Rscript
# Aggregate the archived 1-round vs 2-round effect-modification sweep and plot
# arm-wise TATE bias/RMSE/coverage/CI length.  The original treated-mean figure
# is preserved separately as sim_2round_em_cont_mu1_legacy.pdf.
suppressPackageStartupMessages(library(ggplot2))
input_directory <- "diagnosis/face_probe/validation/chunks2r_cont"
files <- list.files(
  input_directory, pattern = "^em_K4_C1_.*[.]csv$", full.names = TRUE
)
if (!length(files)) stop("no continuous-outcome EM chunk CSVs found")
parts <- lapply(files, function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(error) {
      warning("skipping unreadable chunk ", path, ": ", conditionMessage(error))
      NULL
    }
  )
})
parts <- Filter(Negate(is.null), parts)
if (!length(parts)) stop("no readable continuous-outcome EM chunks found")
raw <- do.call(rbind, parts)
keep <- c(
  "one_round_crossfit_ate", "two_round_crossfit_ate", "target_only_ate"
)
required_columns <- c(
  "sim_id", "method", "em", "p", "bias", "coverage", "ci_width"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("EM chunks are missing columns: ", paste(missing_columns, collapse = ", "))
}
raw <- raw[raw$method %in% keep, required_columns]
if (!nrow(raw)) stop("EM chunks contain no arm-wise TATE rows")
replicate_key <- interaction(
  raw[c("p", "em", "method", "sim_id")],
  drop = TRUE, lex.order = TRUE
)
if (anyDuplicated(replicate_key)) {
  stop("EM chunks contain duplicated p/em/method/sim_id TATE rows")
}

means <- aggregate(
  cbind(bias, coverage, ci_width) ~ method + em + p, raw, mean
)
rmse <- aggregate(
  bias ~ method + em + p, raw, function(value) sqrt(mean(value^2))
)
names(rmse)[names(rmse) == "bias"] <- "rmse"
counts <- aggregate(bias ~ method + em + p, raw, length)
names(counts)[names(counts) == "bias"] <- "n"
agg <- merge(means, rmse, by = c("method", "em", "p"), sort = FALSE)
agg <- merge(agg, counts, by = c("method", "em", "p"), sort = FALSE)
lab <- c(
  one_round_crossfit_ate = "1-round",
  two_round_crossfit_ate = "2-round",
  target_only_ate = "Target-only"
)
agg$Method <- factor(lab[agg$method], levels=c("Target-only","1-round","2-round"))
agg$Pf <- factor(paste0("p = ", agg$p), levels=c("p = 10","p = 50"))
mk <- function(metric, value) data.frame(em=agg$em, Pf=agg$Pf, Method=agg$Method, metric=metric, value=value)
long <- rbind(mk("Bias",agg$bias), mk("RMSE",agg$rmse),
              mk("Coverage",agg$coverage), mk("CI length",agg$ci_width))
long$metric <- factor(
  long$metric, levels = c("Bias", "RMSE", "Coverage", "CI length")
)
href <- data.frame(metric=factor(c("Bias","Coverage"),levels=levels(long$metric)), yint=c(0,0.95))
col_map <- c("Target-only"="#1B9E77","1-round"="#E7298A","2-round"="#D95F02")
lt_map  <- c("Target-only"="longdash","1-round"="dashed","2-round"="solid")
p <- ggplot(long, aes(em, value, colour=Method, linetype=Method, group=Method)) +
  geom_hline(data=href, aes(yintercept=yint), linetype="dashed", colour="grey50", linewidth=0.3) +
  geom_line(linewidth=0.9) +
  facet_grid(metric ~ Pf, scales="free_y") +
  scale_colour_manual(values=col_map) + scale_linetype_manual(values=lt_map) +
  guides(colour=guide_legend(nrow=1), linetype=guide_legend(nrow=1)) +
  labs(x="source effect-modification strength  s", y=NULL) +
  theme_bw(base_size=20) +
  theme(legend.position="top", legend.text=element_text(size=18),
        axis.title=element_text(size=20), axis.text=element_text(size=17),
        strip.text=element_text(size=19),
        panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"))
output_paths <- c(
  "diagnosis/face_probe/validation/fig_2round_em_cont.pdf",
  "docs/figures/sim_2round_em_cont.pdf",
  "overleaf/figures/sim_2round_em_cont.pdf"
)
master_path <- tempfile(
  pattern = ".fig_2round_em_cont_master_",
  tmpdir = dirname(output_paths[[1L]]), fileext = ".pdf"
)
ggsave(master_path, p, width = 8, height = 8.8)
for (output_path in output_paths) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(output_path), "_"),
    tmpdir = dirname(output_path), fileext = ".pdf"
  )
  if (!file.copy(master_path, temporary_path, overwrite = FALSE)) {
    unlink(c(master_path, temporary_path))
    stop("failed to stage ", output_path)
  }
  if (!file.rename(temporary_path, output_path)) {
    unlink(c(master_path, temporary_path))
    stop("failed to atomically install ", output_path)
  }
}
unlink(master_path)
cat(sprintf("wrote fig_2round_em_cont.pdf (%d chunk files)\n", length(files)))
print(agg[order(agg$p, agg$method, agg$em), c("p","method","em","n","bias","rmse","coverage","ci_width")], row.names=FALSE)

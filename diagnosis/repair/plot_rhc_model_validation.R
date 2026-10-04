#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
review <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(file.exists(file.path(review, "CHECKS_PASSED")),
  startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
rows <- read.csv(file.path(review, "replicates.csv"), stringsAsFactors = FALSE)
summary <- read.csv(file.path(review, "summary.csv"), stringsAsFactors = FALSE)
rows <- rows[rows$method %in% c("RoCE", "Target-only"), ]
records <- list()
for (scenario in unique(rows$scenario)) for (method in unique(rows$label)) {
  selected <- rows[rows$scenario == scenario & rows$label == method, ]
  stopifnot(nrow(selected) == 200L, identical(sort(selected$repeat_id), 1:200))
  saved <- summary[summary$scenario == scenario & summary$method == method, ]
  stopifnot(nrow(saved) == 1L)
  error <- selected$estimate - selected$truth
  rmse <- sqrt(mean(error^2)); rmse_mcse <- sd(error^2) / (2 * rmse * sqrt(length(error)))
  stopifnot(abs(rmse - saved$rmse) < 1e-12, abs(mean(error) - saved$bias) < 1e-12,
    abs(mean(selected$covered) - saved$coverage) < 1e-12)
  radius <- if (method == "Target-only") NA_real_ else unique(selected$source_radius)
  records[[length(records) + 1L]] <- data.frame(scenario = scenario, method = method, radius = radius,
    metric = c("Bias (pp)", "RMSE (pp)", "Coverage"),
    value = c(100 * saved$bias, 100 * rmse, saved$coverage),
    lower = c(100 * (saved$bias - qnorm(.975) * saved$bias_mcse),
              100 * max(0, rmse - qnorm(.975) * rmse_mcse), saved$coverage_lower),
    upper = c(100 * (saved$bias + qnorm(.975) * saved$bias_mcse),
              100 * (rmse + qnorm(.975) * rmse_mcse), saved$coverage_upper))
}
values <- do.call(rbind, records)
values$metric <- factor(values$metric, levels = c("Bias (pp)", "RMSE (pp)", "Coverage"))
labels <- c(O_case_mix = "OR correct / case mix", W_overlap = "Weight correct / overlap",
            W_tail = "Weight correct / tail", W_departure = "One source departure",
            O_moderate_mix = "OR correct / moderate shift")
stopifnot(all(unique(values$scenario) %in% names(labels)))
labels <- labels[names(labels) %in% unique(values$scenario)]
values$scenario <- factor(values$scenario, levels = names(labels), labels = unname(labels))
roce <- values[values$method != "Target-only", ]
target <- values[values$method == "Target-only", ]
reference <- data.frame(metric = factor(c("Bias (pp)", "Coverage"), levels = levels(values$metric)),
                        reference = c(0, .95))
suppressPackageStartupMessages(library(ggplot2))
figure <- ggplot(roce, aes(radius, value)) +
  geom_hline(data = reference, aes(yintercept = reference), colour = "grey45", linetype = "dashed") +
  geom_rect(data = target, aes(xmin = 1.8, xmax = 5.2, ymin = lower, ymax = upper),
    inherit.aes = FALSE, fill = "#1B9E77", alpha = .08) +
  geom_hline(data = target, aes(yintercept = value, colour = "Target-only"), linewidth = .7) +
  geom_line(aes(colour = "RoCE"), linewidth = .75) +
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = "RoCE"), width = .13, linewidth = .6) +
  geom_point(aes(colour = "RoCE"), size = 2.5) +
  facet_grid(metric ~ scenario, scales = "free_y") +
  scale_x_continuous(breaks = c(2, 3, 5), limits = c(1.8, 5.2)) +
  scale_colour_manual(values = c("Target-only" = "#1B9E77", "RoCE" = "#D95F02"),
    breaks = c("Target-only", "RoCE"), name = NULL) +
  labs(x = "Source log-weight truncation radius", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold"), legend.position = "bottom")
dir.create(output, recursive = TRUE)
write.csv(values, file.path(output, "plot_data.csv"), row.names = FALSE)
ggsave(file.path(output, "rhc_model_mc200.pdf"), figure, width = max(4.8, 3.2 * length(labels)), height = 7,
  device = grDevices::pdf, useDingbats = FALSE)
writeLines(c(sprintf("%d fixed RHC-feature population laws,200 repeats each,IDs1–200.", length(labels)),
  "Error bars show Monte Carlo95% intervals: normal for bias,delta approximation for RMSE,Wilson for coverage.",
  "Green line/band: calibrated Target-only. Orange: RoCE at each prescribed radius.",
  "When included,the severe case-mix law is a calibration-moment stress test. The moderate OR law has separately checked population moment limits."),
  file.path(output, "caption.txt"))
writeLines("All plotted metrics reproduce the complete200-repeat results.", file.path(output, "CHECKS_PASSED"))

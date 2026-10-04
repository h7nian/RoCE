#!/usr/bin/env Rscript
# Publication plots of the complete saved panel; no refitting or cell selection.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
data_root <- normalizePath(arguments[1L], mustWork = TRUE)
output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
values <- read.csv(file.path(data_root, "simulation_joint_tate_mc200.csv"), stringsAsFactors = FALSE)
display_scenarios <- c(C1 = "Scenario 1", C2 = "Scenario 2", C3 = "Scenario 3")
display_methods <- c(target_anchor_ate = "Target-only", sample_size_ate = "SS",
  inverse_variance_ate = "IVW", federated_dr_ate = "Federated-DR", pooled_dr_ate = "Pooled-DR",
  one_round_crossfit_ate_joint_tate = "RoCE")
dimension <- unique(values$p)
stopifnot(nrow(values) == 420L, all(values$repeats == 200L),
  setequal(unique(values$method), c(names(display_methods), "target_only_ate")),
  length(dimension) == 1L, dimension %in% c(100L, 200L),
  all(is.finite(values$rmse)), all(values$coverage >= 0 & values$coverage <= 1))
selected <- values[values$method %in% names(display_methods), ]
key <- function(data) paste(data$config, data$K, data$rho, sep = ":")
reference <- selected[selected$method == "target_anchor_ate", ]
stopifnot(nrow(reference) == 60L, !anyDuplicated(key(reference)))
selected$rmse_ratio_to_calibrated_target <- selected$rmse/reference$rmse[match(key(selected), key(reference))]
# These archived columns use ordinary AIPW as reference, not the displayed anchor.
selected[c("rmse_ratio_to_target", "paired_mse_difference", "paired_mse_difference_mcse")] <- NULL
selected$method_label <- unname(display_methods[selected$method])
stopifnot(nrow(selected) == 360L, all(is.finite(selected$rmse_ratio_to_calibrated_target)))
for (scenario in names(display_scenarios)) for (sources in c(2, 4, 6, 8)) for (method in names(display_methods)) {
  rows <- selected[selected$config == scenario & selected$K == sources & selected$method == method, ]
  stopifnot(nrow(rows) == 5L, identical(sort(as.numeric(rows$rho)), c(0, .5, 1, 1.5, 2)))
}
dir.create(output, recursive = TRUE)
write.csv(selected, file.path(output, "simulation_joint_tate_mc200.csv"), row.names = FALSE)
provenance <- jsonlite::fromJSON(file.path(data_root, "simulation_joint_tate_provenance.json"))
provenance$methods <- as.list(display_methods)
provenance$scenario_labels <- as.list(display_scenarios)
provenance$metric_rows <- nrow(selected)
provenance$target_only_definition <- "Calibrated target anchor (target_anchor_ate)"
provenance$source_csv_sha256 <- digest::digest(file.path(data_root, "simulation_joint_tate_mc200.csv"), file = TRUE, algo = "sha256")
provenance$display_scope <- "All prespecified cells and replicates; ordinary target-only archived, not displayed"
provenance$figure_style <- "Shared ENAR publication style: Okabe-Ito colours, method shapes and line types"
jsonlite::write_json(provenance, file.path(output, "simulation_joint_tate_provenance.json"), pretty = TRUE, auto_unbox = TRUE)

suppressPackageStartupMessages(library(ggplot2))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(script_path)), "publication_plot_style.R"))
stopifnot(identical(unname(display_methods), publication_method_order))
selected$method_label <- factor(selected$method_label, levels = unname(display_methods))
selected$config <- factor(selected$config, levels = names(display_scenarios))
selected$source_count <- factor(paste0("K = ", selected$K), levels = paste0("K = ", c(2, 4, 6, 8)))
for (metric in c("rmse", "coverage")) {
  selected$value <- selected[[metric]]
  selected$lower <- selected[[paste0(metric, "_lower")]]
  selected$upper <- selected[[paste0(metric, "_upper")]]
  figure <- ggplot(selected, aes(rho, value, colour = method_label, linetype = method_label,
    shape = method_label, group = method_label))
  if (metric == "coverage") {
    figure <- figure + geom_hline(yintercept = .95, colour = "#7D858D", linetype = "dashed", linewidth = .35)
  }
  figure <- figure +
    geom_line(linewidth = .55) +
    geom_point(size = 1.65, stroke = .4) +
    geom_errorbar(data = selected[selected$method == "one_round_crossfit_ate_joint_tate", ],
      aes(ymin = lower, ymax = upper), width = .07, linewidth = .4, linetype = "solid") +
    facet_grid(source_count ~ config, labeller = labeller(config = display_scenarios)) +
    scale_colour_manual(values = publication_colors, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = publication_line_types, name = NULL, drop = FALSE) +
    scale_shape_manual(values = publication_shapes, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(breaks = c(0, .5, 1, 1.5, 2), labels = c("0", "0.5", "1", "1.5", "2")) +
    labs(x = expression("Source deviation " * rho * " (log-odds)"),
      y = if (metric == "rmse") "RMSE" else "Coverage") +
    theme_publication()
  if (metric == "coverage") {
    figure <- figure + scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1))
  } else {
    figure <- figure + scale_y_continuous(limits = c(0, NA))
  }
  ggsave(file.path(output, paste0(metric, "_joint_tate_bounded_v3_mc200.pdf")), figure,
    width = 6.5, height = 6.6, device = grDevices::pdf, useDingbats = FALSE)
}
writeLines(c("ENAR_PUBLICATION_SIMULATION_FIGURES_PASSED",
  "360 displayed metric rows per dimension; all200 replicates in each cell retained.",
  "Target-only is calibrated target_anchor_ate. Orange bars are RoCE Monte Carlo intervals."),
  file.path(output, "checks.txt"))

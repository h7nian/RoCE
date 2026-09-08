# ============================================================================
# plotting.R — ggplot2-based figures for RoCE experiments
# ============================================================================
# Public helpers used by both the simulation pipeline (main.R) and the real-
# data pipeline (realdata.R) to produce publication-ready figures from the
# summary CSVs / RDS files those pipelines emit.
#
# ggplot2 and scales are declared under Suggests in DESCRIPTION so that the
# core package stays lean for environments where graphics are unnecessary
# (e.g., SLURM batch simulation runs on headless compute nodes). Every
# plotting function short-circuits with an informative error when the
# optional packages are absent.
# ============================================================================

utils::globalVariables(".data")

# Internal: guard gracefully when ggplot2 is not installed.
.require_ggplot2 <- function(fun_name = "plotting") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(sprintf(
      "Package 'ggplot2' is required for %s. Install it via install.packages('ggplot2').",
      fun_name
    ))
  }
  invisible(TRUE)
}


# Internal: provide a consistent RoCE colour theme across figures. Uses a
# colour-blind-safe palette (Okabe-Ito) and keeps RoCE highlighted.
.roce_colours <- function(methods, highlight = "roce") {
  # Reserve vermillion for the highlighted method. Keeping it out of the base
  # cycle prevents a neighbouring method from receiving the same colour when
  # factor levels are reversed in a forest plot.
  palette <- c(
    "#0072B2", "#009E73", "#CC79A7", "#56B4E9",
    "#E69F00", "#F0E442", "#999999", "#000000"
  )
  out <- setNames(rep(palette, length.out = length(methods)), methods)
  if (!is.null(highlight) && highlight %in% methods) {
    out[highlight] <- "#D55E00"  # vermillion accent for RoCE
  }
  out
}


# Internal: default theme shared by all RoCE figures.  The manuscript uses
# multi-panel, full-width figures, so 20 pt is the minimum readable default
# after LaTeX scaling. Callers can still override it for a different layout.
.roce_theme <- function(base_size = 20) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey92"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      strip.text = ggplot2::element_text(
        face = "bold", size = base_size
      ),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        face = "bold", size = base_size - 1
      ),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key = ggplot2::element_blank()
    )
}


# ============================================================================
# Shared helper: save a ggplot to disk.
# ============================================================================

#' Save a RoCE ggplot to PDF / PNG
#'
#' Thin wrapper around \code{ggplot2::ggsave} that uses RoCE defaults for
#' width, height, and dpi. Returns the output path invisibly so call sites can
#' chain pipes.
#'
#' @param plot A \code{ggplot} object.
#' @param path Output file path. Must end in \code{.pdf} or \code{.png}.
#' @param width,height Numeric (inches). Defaults tuned for two-column
#'   journal figures (7 x 4.5 inches).
#' @param dpi Integer. Defaults to 300; applied only to PNG output.
#' @return \code{path}, returned invisibly.
#' @export
save_plot <- function(plot, path, width = 7, height = 4.5, dpi = 300) {
  .require_ggplot2("save_plot()")
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.")
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = path, plot = plot,
                  width = width, height = height, dpi = dpi, units = "in")
  invisible(path)
}


# ============================================================================
# Real-data plots
# ============================================================================

#' Forest Plot of Method-wise Point Estimates and 95\% CIs
#'
#' Horizontal point-and-interval plot comparing the target-effect estimate of
#' several methods (target-only, SS, IVW, Tilted-AIPW, RoCE, ...). Suitable
#' for the real-data application section. Accepts either a long data frame
#' with the standard RoCE comparison-method columns, or the list returned by
#' \code{run_rhc_experiment()}.
#'
#' @param methods_df Data frame with required columns \code{method},
#'   \code{estimate}, \code{ci_lower}, \code{ci_upper}. A \code{se} column is
#'   used when CIs are missing.
#' @param highlight Character. Method label to visually accent (default
#'   \code{"roce"}).
#' @param title,subtitle Optional plot annotations.
#' @param xlab Custom x-axis label; default \code{"Treatment effect estimate"}.
#' @return A \code{ggplot} object.
#' @export
plot_forest_methods <- function(methods_df,
                                highlight = "roce",
                                title = NULL,
                                subtitle = NULL,
                                xlab = "Treatment effect estimate") {
  .require_ggplot2("plot_forest_methods()")
  required <- c("method", "estimate")
  missing_cols <- setdiff(required, colnames(methods_df))
  if (length(missing_cols) > 0L) {
    stop("methods_df is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  df <- as.data.frame(methods_df)
  if (!all(c("ci_lower", "ci_upper") %in% colnames(df))) {
    if (!"se" %in% colnames(df)) {
      stop("Either (ci_lower, ci_upper) or se must be provided in methods_df.")
    }
    df$ci_lower <- df$estimate - 1.96 * df$se
    df$ci_upper <- df$estimate + 1.96 * df$se
  }

  df$method <- factor(df$method, levels = rev(unique(df$method)))
  colours <- .roce_colours(levels(df$method), highlight = highlight)

  ggplot2::ggplot(df, ggplot2::aes(y = .data$method, x = .data$estimate)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper,
                   colour = .data$method),
      height = 0.25, linewidth = 0.8
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$method),
                        size = 2.8) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    .roce_theme()
}


#' Bar Chart of RoCE Aggregation Weights
#'
#' Displays the learned common TATE source weights
#' \eqn{\widehat{\eta}_j} returned by \code{\link{run_tate_crossfit}}. A zero
#' weight is an optimizer outcome; crossing the Wald activation cutoff alone
#' does not guarantee an exact zero in finite samples.
#'
#' @param weights Numeric vector of length \eqn{K} (source-site weights).
#' @param source_labels Optional character vector of labels (default
#'   \code{s1, s2, ...}).
#' @param penalty_d2 Deprecated optional squared-discrepancy annotations retained
#'   only for compatibility with archived treated-mean plotting code. Supply at
#'   most one of \code{wald_statistics} and \code{penalty_d2}.
#' @param title,subtitle Optional annotations.
#' @param wald_statistics Optional numeric vector of source-specific TATE
#'   Wald discrepancy statistics, used to annotate the bars. This argument is
#'   placed last to preserve the positional API of earlier releases.
#' @return A \code{ggplot} object.
#' @export
plot_aggregation_weights <- function(weights,
                                     source_labels = NULL,
                                     penalty_d2 = NULL,
                                     title = NULL,
                                     subtitle = NULL,
                                     wald_statistics = NULL) {
  .require_ggplot2("plot_aggregation_weights()")
  if (!is.numeric(weights)) stop("weights must be numeric.")
  K <- length(weights)
  if (is.null(source_labels)) source_labels <- paste0("s", seq_len(K))
  if (length(source_labels) != K) {
    stop("length(source_labels) must equal length(weights).")
  }
  if (!is.null(wald_statistics) && !is.null(penalty_d2)) {
    stop("Supply at most one of wald_statistics and penalty_d2.")
  }

  df <- data.frame(
    source    = factor(source_labels, levels = source_labels),
    weight    = weights,
    nonzero   = abs(weights) > .Machine$double.eps,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(
      df, ggplot2::aes(x = .data$source, y = .data$weight,
                       fill = .data$nonzero)
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#0072B2", `FALSE` = "grey70"),
      labels = c(`TRUE` = "Nonzero", `FALSE` = "Zero"),
      name = "Weight status"
    ) +
    ggplot2::labs(
      x = "Source site", y = expression(hat(eta)[j]),
      title = title, subtitle = subtitle
    ) +
    .roce_theme()

  if (!is.null(wald_statistics)) {
    if (!is.numeric(wald_statistics) || length(wald_statistics) != K ||
        any(!is.finite(wald_statistics)) || any(wald_statistics < 0)) {
      stop("wald_statistics must contain K finite non-negative values.")
    }
    df$wald_statistic <- wald_statistics
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(label = sprintf("t==%.2f", .data$wald_statistic)),
      parse = TRUE, vjust = -0.5, size = 5
    )
  } else if (!is.null(penalty_d2)) {
    if (length(penalty_d2) != K) stop("penalty_d2 must have length K.")
    df$d2 <- penalty_d2
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(label = sprintf("d^2==%.3g", .data$d2)),
      parse = TRUE, vjust = -0.5, size = 5
    )
  }
  p
}


#' Scatter of Pairwise vs. Target-Only Estimates
#'
#' Visualizes each source-assisted pairwise estimate
#' \eqn{\widehat{\mu}^a_{t,s_j}} against the target-only estimate
#' \eqn{\widehat{\mu}^a_{ot}}. Points close to the \eqn{y = x} reference line
#' indicate informative sources that the RoCE aggregation is expected to
#' retain; points far from the line correspond to sources that the selection
#' penalty should drive toward zero.
#'
#' @param target_only_estimate Numeric scalar, \eqn{\widehat{\mu}^a_{ot}}.
#' @param pairwise Data frame with columns \code{source},
#'   \code{estimate}, and optionally \code{se} (used for error bars).
#' @param title,subtitle Optional annotations.
#' @return A \code{ggplot} object.
#' @export
plot_pairwise_vs_target <- function(target_only_estimate,
                                    pairwise,
                                    title = NULL,
                                    subtitle = NULL) {
  .require_ggplot2("plot_pairwise_vs_target()")
  required <- c("source", "estimate")
  if (!all(required %in% colnames(pairwise))) {
    stop("pairwise must contain columns ", paste(required, collapse = ", "))
  }

  df <- as.data.frame(pairwise)
  df$target <- target_only_estimate
  has_se <- "se" %in% colnames(df)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$target, y = .data$estimate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(size = 2.8, colour = "#0072B2") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$source),
      hjust = -0.3, vjust = 0, size = 3
    ) +
    ggplot2::labs(
      x = expression(hat(mu)[ot]^a),
      y = expression(hat(mu)["t,s"[j]]^a),
      title = title, subtitle = subtitle
    ) +
    .roce_theme()

  if (has_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$estimate - 1.96 * .data$se,
                   ymax = .data$estimate + 1.96 * .data$se),
      width = 0.0, colour = "#0072B2", alpha = 0.5
    )
  }
  p
}


# ============================================================================
# Simulation plots
# ============================================================================

#' Method-comparison Plot for Simulation Summary Output
#'
#' Produces a metric-faceted bar chart for the output of
#' \code{\link{summarize_results}} across simulation configurations. The
#' default metric set (\code{bias_mean}, \code{bias_rmse}, \code{coverage},
#' \code{ci_width}) matches what the simulation pipeline writes to
#' \code{results/}.
#'
#' @param summary_df Data frame from \code{\link{summarize_results}}, read
#'   back from the CSV written by \code{main.R}. Expected columns:
#'   \code{method}, \code{config}, \code{n_total}, \code{K}, and one column
#'   per metric (e.g. \code{bias.mean}, \code{bias.rmse}, \code{coverage},
#'   \code{ci_width}).
#' @param metric One of \code{"bias_mean"}, \code{"bias_rmse"},
#'   \code{"coverage"}, \code{"ci_width"}.
#' @param facet_cols Character vector of columns on which to facet. Defaults
#'   to \code{c("config")}.
#' @param x Character. Variable used on the x-axis. Defaults to \code{"n_total"}.
#' @param highlight Method label to accent (default \code{"roce"}).
#' @param title,subtitle Optional annotations.
#' @return A \code{ggplot} object.
#' @export
plot_simulation_metric <- function(summary_df,
                                   metric = c("bias_mean", "bias_rmse",
                                              "coverage", "ci_width"),
                                   facet_cols = "config",
                                   x = "n_total",
                                   highlight = "roce",
                                   title = NULL,
                                   subtitle = NULL) {
  .require_ggplot2("plot_simulation_metric()")
  metric <- match.arg(metric)

  # Normalise metric-column names across the two CSV conventions used in
  # main.R (bias.mean vs bias_mean, depending on aggregate() output).
  metric_col <- switch(
    metric,
    bias_mean = if ("bias.mean" %in% colnames(summary_df)) "bias.mean" else "bias_mean",
    bias_rmse = if ("bias.rmse" %in% colnames(summary_df)) "bias.rmse" else "bias_rmse",
    coverage  = "coverage",
    ci_width  = "ci_width"
  )
  if (!metric_col %in% colnames(summary_df)) {
    stop("metric column not found in summary_df: ", metric_col)
  }

  df <- as.data.frame(summary_df)
  df$method <- as.character(df$method)
  methods   <- unique(df$method)
  colours   <- .roce_colours(methods, highlight = highlight)

  y_label <- switch(
    metric,
    bias_mean = "Mean bias",
    bias_rmse = "RMSE",
    coverage  = "95% CI coverage",
    ci_width  = "Mean CI width"
  )

  # Use tidy-evaluation `.data[[col]]` lookups instead of the deprecated
  # aes_string(); this keeps string-valued column names (from the CSV) while
  # remaining compatible with ggplot2 >= 3.0 without deprecation warnings.
  p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[x]], y = .data[[metric_col]],
                   colour = .data[["method"]], group = .data[["method"]])
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = colours, name = "Method") +
    ggplot2::labs(x = x, y = y_label, title = title, subtitle = subtitle) +
    .roce_theme()

  # Add a 0.95 reference line for coverage plots.
  if (metric == "coverage") {
    p <- p + ggplot2::geom_hline(yintercept = 0.95,
                                 linetype = "dashed", colour = "grey40")
  }

  if (length(facet_cols) > 0L) {
    facet_formula <- stats::as.formula(paste("~", paste(facet_cols, collapse = " + ")))
    p <- p + ggplot2::facet_wrap(facet_formula, scales = "free_y")
  }
  p
}


#' Convenience: produce and save the full simulation plot grid
#'
#' Reads a summary CSV written by \code{main.R}, builds all four metric
#' plots (\code{bias_mean}, \code{bias_rmse}, \code{coverage},
#' \code{ci_width}), and saves them as PDFs under \code{out_dir}.
#'
#' @param summary_csv Path to the summary CSV produced by \code{main.R}.
#' @param out_dir Directory to save the four PDFs. Created if absent.
#' @param facet_cols,x,highlight Passed to \code{\link{plot_simulation_metric}}.
#' @return Invisible character vector of written file paths.
#' @export
save_simulation_plot_grid <- function(summary_csv,
                                      out_dir,
                                      facet_cols = "config",
                                      x = "n_total",
                                      highlight = "roce") {
  .require_ggplot2("save_simulation_plot_grid()")
  if (!file.exists(summary_csv)) {
    stop("summary_csv not found: ", summary_csv)
  }
  summary_df <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)

  metrics <- c("bias_mean", "bias_rmse", "coverage", "ci_width")
  paths <- character(length(metrics))
  for (i in seq_along(metrics)) {
    p <- plot_simulation_metric(summary_df, metric = metrics[i],
                                facet_cols = facet_cols, x = x,
                                highlight = highlight)
    paths[i] <- save_plot(
      p,
      file.path(out_dir, sprintf("simulation_%s.pdf", metrics[i])),
      width = 8, height = 5
    )
  }
  invisible(paths)
}

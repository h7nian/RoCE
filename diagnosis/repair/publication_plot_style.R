# Shared presentation for the ENAR simulation and application figures.
# Colours follow Okabe--Ito; shapes and line types also identify the methods.
publication_method_order <- c("Target-only", "SS", "IVW", "Federated-DR", "Pooled-DR", "RoCE")
publication_colors <- c("Target-only" = "#444444", "SS" = "#E69F00", "IVW" = "#56B4E9",
  "Federated-DR" = "#009E73", "Pooled-DR" = "#CC79A7", "RoCE" = "#D55E00")
publication_line_types <- c("Target-only" = "longdash", "SS" = "dotted", "IVW" = "dashed",
  "Federated-DR" = "dotdash", "Pooled-DR" = "twodash", "RoCE" = "solid")
publication_shapes <- c("Target-only" = 16, "SS" = 17, "IVW" = 15,
  "Federated-DR" = 18, "Pooled-DR" = 8, "RoCE" = 19)

theme_publication <- function(base_size = 10.5) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#292929"),
      axis.text = ggplot2::element_text(colour = "#444444"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 7)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E4E7EB", linewidth = .3),
      panel.border = ggplot2::element_rect(colour = "#D6DADF", fill = NA, linewidth = .3),
      panel.spacing = grid::unit(9, "pt"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(21, "pt"),
      legend.key.height = grid::unit(14, "pt"),
      legend.margin = ggplot2::margin(0, 0, 5, 0),
      plot.margin = ggplot2::margin(6, 8, 6, 6))
}

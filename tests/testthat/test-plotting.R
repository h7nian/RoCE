# Tests for R/plotting.R. The package declares ggplot2 under Suggests, so
# every test is skipped gracefully when the optional dependency is absent.

skip_if_no_ggplot2 <- function() {
  skip_if_not_installed("ggplot2")
}

.fake_methods_df <- function() {
  data.frame(
    method   = c("target_only", "sample_size", "inverse_variance",
                 "tilted_aipw", "facehd"),
    estimate = c(0.60, 0.55, 0.58, 0.62, 0.59),
    se       = c(0.04, 0.03, 0.03, 0.04, 0.02),
    ci_lower = c(0.52, 0.49, 0.52, 0.54, 0.55),
    ci_upper = c(0.68, 0.61, 0.64, 0.70, 0.63),
    stringsAsFactors = FALSE
  )
}

test_that("plot_forest_methods returns a ggplot object", {
  skip_if_no_ggplot2()
  p <- plot_forest_methods(.fake_methods_df(), highlight = "facehd")
  expect_s3_class(p, "ggplot")
})

test_that("plot_forest_methods derives CIs from se when missing", {
  skip_if_no_ggplot2()
  df <- .fake_methods_df()[, c("method", "estimate", "se")]  # drop CIs
  p <- plot_forest_methods(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_forest_methods errors on missing required columns", {
  skip_if_no_ggplot2()
  expect_error(plot_forest_methods(data.frame(method = "facehd")),
               "missing required columns")
})

test_that("plot_aggregation_weights returns a ggplot object", {
  skip_if_no_ggplot2()
  w <- c(0.3, 0.0, 0.25, 0.0)
  p <- plot_aggregation_weights(w, source_labels = paste0("s", 1:4),
                                penalty_d2 = c(0.01, 0.04, 0.005, 0.09))
  expect_s3_class(p, "ggplot")
})

test_that("plot_pairwise_vs_target returns a ggplot object", {
  skip_if_no_ggplot2()
  pairwise <- data.frame(
    source   = paste0("s", 1:4),
    estimate = c(0.55, 0.62, 0.49, 0.71),
    se       = c(0.04, 0.03, 0.05, 0.02),
    stringsAsFactors = FALSE
  )
  p <- plot_pairwise_vs_target(0.60, pairwise)
  expect_s3_class(p, "ggplot")
})

test_that("plot_simulation_metric handles both bias.mean and bias_mean columns", {
  skip_if_no_ggplot2()
  df_dot <- data.frame(
    method    = rep(c("facehd", "target_only"), each = 3),
    config    = "C1",
    n_total   = rep(c(5000, 10000, 20000), 2),
    K         = 3,
    bias.mean = c(0.01, 0.008, 0.005, 0.04, 0.03, 0.02),
    stringsAsFactors = FALSE
  )
  expect_s3_class(
    plot_simulation_metric(df_dot, metric = "bias_mean"),
    "ggplot"
  )

  df_under <- df_dot
  colnames(df_under)[colnames(df_under) == "bias.mean"] <- "bias_mean"
  expect_s3_class(
    plot_simulation_metric(df_under, metric = "bias_mean"),
    "ggplot"
  )
})

test_that("save_plot writes a PDF and returns the path invisibly", {
  skip_if_no_ggplot2()
  tmp <- tempfile(fileext = ".pdf")
  p <- plot_forest_methods(.fake_methods_df())
  out <- save_plot(p, tmp, width = 5, height = 3)
  expect_true(file.exists(tmp))
  expect_identical(out, tmp)
  unlink(tmp)
})

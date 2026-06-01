# Tests for the Right-Heart Catheterization (RHC) public-dataset loader.
# Skipped automatically when inst/extdata/rhc.csv has not yet been vendored
# (e.g., a fresh clone before running data/download_rhc.sh).

skip_if_no_rhc_csv <- function() {
  path <- rhc_csv_path()
  if (!nzchar(path) || !file.exists(path)) {
    skip("inst/extdata/rhc.csv not vendored; run data/download_rhc.sh first.")
  }
}

test_that("rhc_csv_path returns a string", {
  expect_type(rhc_csv_path(), "character")
})

test_that("load_rhc_raw returns a non-empty data frame with expected columns", {
  skip_if_no_rhc_csv()
  raw <- load_rhc_raw()
  expect_s3_class(raw, "data.frame")
  expect_gt(nrow(raw), 5000L)
  required_cols <- c("cat1", "swang1", "dth30", "death", "age", "sex",
                     "aps1", "scoma1", "meanbp1", "wblc1", "ninsclas")
  expect_true(all(required_cols %in% colnames(raw)))
})

test_that("build_rhc_cohort returns a coherent frame with A and Y", {
  skip_if_no_rhc_csv()
  cohort <- build_rhc_cohort(outcome = "death30")

  expect_s3_class(cohort, "data.frame")
  expect_true(all(c("A", "Y", "site_var") %in% colnames(cohort)))
  expect_true(all(cohort$A %in% c(0L, 1L)))
  expect_true(all(cohort$Y %in% c(0L, 1L)))
  expect_equal(attr(cohort, "outcome"), "death30")

  # Covariate columns are finite (no NAs left after imputation).
  covariates <- setdiff(colnames(cohort), c("A", "Y", "site_var"))
  expect_true(length(covariates) > 40L)
  for (col in covariates) {
    expect_true(all(is.finite(cohort[[col]])),
                info = sprintf("Column '%s' contains non-finite values", col))
  }
})

test_that("build_rhc_cohort supports binary death180 and continuous los outcomes", {
  skip_if_no_rhc_csv()
  cohort_bin <- build_rhc_cohort(outcome = "death180")
  expect_true(all(cohort_bin$Y %in% c(0L, 1L)))
  expect_equal(attr(cohort_bin, "outcome"), "death180")

  cohort_los <- build_rhc_cohort(outcome = "los")
  expect_true(is.numeric(cohort_los$Y))
  expect_true(all(is.finite(cohort_los$Y)))
  expect_equal(attr(cohort_los, "outcome"), "los")
})

test_that("build_rhc_data_split returns a valid FACE-HD data_split", {
  skip_if_no_rhc_csv()
  data_split <- build_rhc_data_split(K = 5L)

  # Naming convention: "t" + "s1", ..., "s<K-1>"
  expect_equal(names(data_split)[1L], "t")
  expect_equal(length(data_split), 5L)
  source_labels <- setdiff(names(data_split), "t")
  expect_equal(source_labels, paste0("s", seq_len(4L)))

  required_fields <- c("n", "X", "X_dagger", "A", "Y", "Z_site", "W_outcome")
  for (site in names(data_split)) {
    site_data <- data_split[[site]]
    expect_true(all(required_fields %in% names(site_data)),
                info = sprintf("site '%s' is missing required fields", site))
    expect_equal(nrow(site_data$X), site_data$n)
    expect_equal(nrow(site_data$X_dagger), site_data$n)
    expect_equal(length(site_data$A), site_data$n)
    expect_equal(length(site_data$Y), site_data$n)
    expect_equal(nrow(site_data$Z_site), site_data$n)
    expect_equal(nrow(site_data$W_outcome), site_data$n)
    expect_true(all(site_data$A %in% c(0L, 1L)))
  }

  # Default phi = identity means W_outcome matches X column-for-column.
  expect_identical(dim(data_split$t$W_outcome), dim(data_split$t$X))
  expect_false(attr(data_split, "phi_applied"))
})

test_that("build_rhc_data_split accepts a custom phi respecting the interface", {
  skip_if_no_rhc_csv()

  double_phi <- function(X) cbind(X, X)
  data_split <- build_rhc_data_split(K = 3L, phi = double_phi)

  expect_true(attr(data_split, "phi_applied"))
  expect_equal(ncol(data_split$t$W_outcome), 2L * ncol(data_split$t$X))
  # X remains the unexpanded raw covariate matrix.
  expect_lt(ncol(data_split$t$X), ncol(data_split$t$W_outcome))
})

test_that("build_rhc_data_split honours target_site override", {
  skip_if_no_rhc_csv()
  default_ds  <- build_rhc_data_split(K = 5L)
  default_map <- attr(default_ds, "site_mapping")

  alt_target  <- unname(default_map["s1"])  # second-largest site
  custom_ds   <- build_rhc_data_split(K = 5L, target_site = alt_target)
  custom_map  <- attr(custom_ds, "site_mapping")

  expect_equal(unname(custom_map["target"]), alt_target)
  expect_false(identical(default_map["target"], custom_map["target"]))
})

test_that("phi_rhc_bspline returns an expanded matrix and preserves row count", {
  skip_if_no_rhc_csv()
  cohort <- build_rhc_cohort()
  # Use a numeric-only subset of the cohort (drop A, Y, and the character
  # site_var column) so the double-coercion below is lossless.
  numeric_cols <- setdiff(colnames(cohort), c("A", "Y", "site_var"))
  X <- as.matrix(cohort[, numeric_cols[seq_len(10L)], drop = FALSE])
  storage.mode(X) <- "double"

  X_expanded <- phi_rhc_bspline(X, knots_per_var = 3L, degree = 3L)

  expect_true(is.matrix(X_expanded))
  expect_equal(nrow(X_expanded), nrow(X))
  expect_gte(ncol(X_expanded), ncol(X))
})

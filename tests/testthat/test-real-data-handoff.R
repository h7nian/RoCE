test_that("real-data adapter preserves labels, rows and shared feature maps", {
  helper <- file.path(repo_root, "scripts/real_data/helpers.R")
  skip_if_not(file.exists(helper), "Repository-only collaborator adapter")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  frame <- data.frame(site = c("other", "target", "third", "target"),
                      A = c(1, 0, 0, 1), Y = c(0, 1, 0, 1), x1 = 1:4, x2 = 4:1)
  build <- function(data = frame, features = c("x2", "x1")) {
    scope$site_data_from_frame(data, "site", "A", "Y", features, "target")
  }
  result <- build()
  expect_identical(names(result), c("t", "s1", "s2"))
  expect_identical(unname(attr(result, "site_mapping")), c("target", "other", "third"))
  expect_identical(result$t$A, frame$A[c(2, 4)])
  expected <- as.matrix(frame[c(2, 4), c("x2", "x1")])
  rownames(expected) <- NULL
  expect_identical(result$t$W_outcome, expected)
  expect_identical(result$t$W_outcome, result$t$Z_site)
  expect_error(build(features = c("x1", "A")), "distinct")
  frame$x1[1] <- NA_real_
  expect_error(build(frame), "finite")
  frame$x1 <- factor(1:4)
  expect_error(build(frame), "numeric")
})

test_that("real-data curvature output contains only aggregate solver diagnostics", {
  scope <- new.env(parent = globalenv())
  helper <- file.path(repo_root, "scripts/real_data/helpers.R")
  skip_if_not(file.exists(helper), "Repository-only collaborator adapter")
  sys.source(helper, envir = scope)
  design <- cbind(1, 1:10, 1:10)
  solved <- .solve_baseline_curvature(design, rep(1, 10), colMeans(design), list(score = design))
  fit <- list(components = list(curvature_diagnostics = list(
    mu1 = list(s1 = list(PS = solved$diagnostics, OR = NULL)), mu0 = NULL)))
  rows <- scope$real_data_curvature_rows(fit)
  expect_equal(nrow(rows), 1L)
  expect_identical(rows$solver, "identifiable_subspace")
  expect_identical(rows$site, "s1")
  expect_identical(rows$arm, "mu1")
  expect_equal(rows$active_columns, 3)
  expect_equal(rows$rank, 2)
  expect_false(any(grepl("patient|coefficient|feature_name|prediction", names(rows))))
})

test_that("real-data runner persists warnings and failure status before exiting", {
  runner <- file.path(repo_root, "scripts/real_data/run.R")
  skip_if_not(file.exists(runner), "Repository-only collaborator runner")
  work <- tempfile("runner_failure_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  set.seed(1169)
  x <- matrix(rnorm(120 * 4), 120, 4, dimnames = list(NULL, paste0("x", 1:4)))
  site <- list(X = x, X_dagger = x, W_outcome = x, Z_site = x,
    A = rep(0:1, 60), Y = rep(c(0, 0, 1, 1), 30), n = 120L)
  prepared <- list(data_split = list(t = site, s1 = site), comparison_data = list(t = site, s1 = site),
    folds = list(target_folds = list()), metadata = list(synthetic = TRUE))
  saveRDS(prepared, file.path(work, "prepared.rds"))
  literal <- function(x) paste(capture.output(dput(x)), collapse = "\n")
  config <- file.path(work, "config.R")
  writeLines(c("list(seed = 11L, n_bootstrap = 100L,",
    "fit = list(family = 'binomial', n_folds = 3L),",
    paste0("build_data = function() readRDS(", literal(file.path(work, "prepared.rds")), "))")), config)
  output <- file.path(work, "output")
  wrapper <- file.path(work, "failure.R")
  # Inject a warning followed by an error at the fit boundary, in an isolated
  # R process. The real runner must save both and must never mark completion.
  writeLines(c(
    paste0(".libPaths(", literal(.libPaths()), ")"),
    "library(RoCE)",
    "trace('.target_tate_reference', where = asNamespace('RoCE'), print = FALSE,",
    "  tracer = quote({ warning('runner warning sentinel'); stop('runner failure sentinel') }))",
    "scope <- new.env(parent = globalenv())",
    paste0("scope$commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) ",
      literal(c(config, output, "federated_dr")), " else ", literal(paste0("--file=", runner))),
    paste0("sys.source(", literal(runner), ", envir = scope)")), wrapper)
  log <- file.path(work, "run.log")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(wrapper)),
    stdout = log, stderr = log, env = "ROCE_CORES=1")
  expect_true(status != 0L)
  directory <- file.path(output, "federated_dr")
  expect_true(file.exists(file.path(directory, "warnings.txt")))
  expect_true(file.exists(file.path(directory, "status.txt")))
  expect_match(paste(readLines(file.path(directory, "warnings.txt")), collapse = "\n"),
    "runner warning sentinel")
  status_text <- paste(readLines(file.path(directory, "status.txt")), collapse = "\n")
  expect_match(status_text, "status: failed")
  expect_match(status_text, "method: federated_dr")
  expect_match(status_text, "runner failure sentinel")
  expect_false(file.exists(file.path(directory, "COMPLETE")))
})

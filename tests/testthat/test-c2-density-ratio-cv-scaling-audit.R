test_that("density-ratio CV validation loss uses full-source empirical scaling", {
  skip_if_not(
    Sys.getenv("FACEHD_C2_CV_SCALE_AUDIT", "0") %in% c("1", "TRUE", "true", "True"),
    "Set FACEHD_C2_CV_SCALE_AUDIT=1 to run the C2 density-ratio CV scaling audit."
  )

  repo_root <- normalizePath(file.path(dirname(test_path()), "..", ".."), mustWork = TRUE)
  cv_utils_path <- file.path(repo_root, "src", "cv_utils.hpp")
  density_path <- file.path(repo_root, "src", "density_ratio.cpp")

  cv_utils <- paste(readLines(cv_utils_path, warn = FALSE), collapse = "\n")
  density_cpp <- paste(readLines(density_path, warn = FALSE), collapse = "\n")

  expect_match(
    cv_utils,
    "density_ratio_cd_update[\\s\\S]*int n_total[\\s\\S]*grad_acc\\.add\\([^\\n]+/ n_total\\)",
    perl = TRUE,
    info = "Training loss should keep the source empirical expectation on the full source-sample scale."
  )

  val_signature <- regmatches(
    cv_utils,
    regexpr("density_ratio_val_loss\\([^\\)]*\\)", cv_utils, perl = TRUE)
  )
  expect_true(length(val_signature) == 1L && nzchar(val_signature))

  expect_match(
    val_signature,
    "source_scale|arm_fraction|n_source|n_total",
    perl = TRUE,
    info = paste(
      "The validation loss is evaluated on treated validation folds.",
      "It needs an explicit full-source scaling factor such as n_treated / n_source."
    )
  )

  expect_false(
    grepl("exp_neg_g * psi_prime_val(j) / n_val", cv_utils, fixed = TRUE),
    info = paste(
      "A treated-fold mean estimates E[f(X) | A=a], not",
      "E[I(A=a) f(X)] on the source empirical scale used by training and main.tex."
    )
  )

  val_calls <- gregexpr("density_ratio_val_loss\\([^;]+\\)", density_cpp, perl = TRUE)
  val_call_text <- regmatches(density_cpp, val_calls)[[1]]
  expect_true(length(val_call_text) >= 3L)
  expect_true(
    all(grepl("source_scale|arm_fraction|n_treated\\s*/\\s*n|n_treated\\s*/\\s*static_cast", val_call_text, perl = TRUE)),
    info = paste(
      "Every density-ratio CV caller should pass the source-arm fraction",
      "to keep validation and training losses on the same empirical scale."
    )
  )
})

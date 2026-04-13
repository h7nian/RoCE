# This file is part of the testthat test framework for FACE-C
# Run all tests with: testthat::test_dir("tests/testthat")

library(testthat)

# Load the FACEC package (prefer installed; fall back to devtools::load_all)
if (requireNamespace("FACEC", quietly = TRUE)) {
  library(FACEC)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  # Determine project root (works from tests/ or project root)
  if (file.exists("../../DESCRIPTION")) {
    devtools::load_all("../..")
  } else if (file.exists("DESCRIPTION")) {
    devtools::load_all(".")
  }
}

# Use test_dir instead of test_check (works without installed package)
test_dir("tests/testthat", reporter = "summary")

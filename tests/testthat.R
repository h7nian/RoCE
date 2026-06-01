# This file is part of the testthat test framework for FACE-HD
# Run all tests with: testthat::test_dir("tests/testthat")

library(testthat)

# Load the FACEHD package for standalone testthat runs.
if (requireNamespace("FACEHD", quietly = TRUE)) {
  library(FACEHD)
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

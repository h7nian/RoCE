# Backward-compatible repository entry point.  The canonical implementation is
# also bundled with package tests so isolated `R CMD check` runs do not depend
# on files excluded by .Rbuildignore.
helper_path <- file.path("tests", "testthat", "helper-c2-true-gamma-utils.R")
if (!file.exists(helper_path)) {
  stop("Cannot locate shared C2 true-gamma helpers at ", helper_path,
       call. = FALSE)
}
source(helper_path, local = FALSE)

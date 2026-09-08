# test-constants_sync.R - Verify R and C++ constants remain in sync
#
# Reads src/numerical_constants.h and compares parsed values against
# the corresponding R constants in R/constants.R.
# Prevents silent drift between the two files.

test_that("R and C++ numerical constants are in sync", {
  # Read and parse C++ constants from header file
  header_path <- file.path("..", "..", "src", "numerical_constants.h")
  if (!file.exists(header_path)) {
    header_path <- file.path(system.file(package = "RoCE"), "..", "src", "numerical_constants.h")
  }
  # Also try the common testthat working directory
  if (!file.exists(header_path)) {
    header_path <- "src/numerical_constants.h"
  }
  skip_if_not(file.exists(header_path), "Cannot locate numerical_constants.h")
  
  header_lines <- readLines(header_path)
  
  # Helper: extract constexpr values from C++ header
  parse_cpp_constant <- function(name) {
    pattern <- sprintf("constexpr\\s+\\w+\\s+%s\\s*=\\s*([^;]+);", name)
    match <- grep(pattern, header_lines, value = TRUE)
    if (length(match) == 0) return(NULL)
    val_str <- sub(paste0(".*", name, "\\s*=\\s*"), "", match[1])
    val_str <- sub(";.*", "", val_str)
    val_str <- trimws(val_str)
    # Handle negative values and scientific notation
    as.numeric(val_str)
  }
  
  # Map of C++ constant name -> R constant name -> expected value
  # These are the constants that MUST be identical across languages
  sync_pairs <- list(
    list(cpp = "M_TAU_DEFAULT",       r_val = M_TAU_DEFAULT),
    list(cpp = "VAR_MIN",             r_val = VARIANCE_MIN),
    list(cpp = "WEIGHT_MIN",          r_val = WEIGHT_MIN),
    list(cpp = "WEIGHT_MAX",          r_val = WEIGHT_MAX),
    list(cpp = "PARAM_MAX",           r_val = PARAM_MAX),
    list(cpp = "CV_FAILURE_PATIENCE", r_val = NUISANCE_CV_FAILURE_PATIENCE),
    list(cpp = "ACTIVE_SET_THRESHOLD", r_val = ACTIVE_SET_THRESHOLD),
    list(cpp = "HESSIAN_FLOOR",       r_val = HESSIAN_FLOOR)
  )
  
  # Also check ETA_CLIP_MIN/MAX vs LOGISTIC_CLIP
  eta_clip_min <- parse_cpp_constant("ETA_CLIP_MIN")
  eta_clip_max <- parse_cpp_constant("ETA_CLIP_MAX")
  if (!is.null(eta_clip_min) && !is.null(eta_clip_max)) {
    expect_equal(abs(eta_clip_min), LOGISTIC_CLIP,
                 info = "ETA_CLIP_MIN magnitude must match R LOGISTIC_CLIP")
    expect_equal(eta_clip_max, LOGISTIC_CLIP,
                 info = "ETA_CLIP_MAX must match R LOGISTIC_CLIP")
  }
  
  for (pair in sync_pairs) {
    cpp_val <- parse_cpp_constant(pair$cpp)
    if (!is.null(cpp_val)) {
      expect_equal(cpp_val, pair$r_val,
                   info = sprintf("C++ %s (%.10g) != R value (%.10g)", 
                                  pair$cpp, cpp_val, pair$r_val))
    }
  }
})

test_that("glmnet path iteration budget is explicit and conservative", {
  expect_identical(GLMNET_MAX_ITER, 1000000L)
  expect_gt(GLMNET_MAX_ITER, 100000L)
})

test_that("production aggregation cutoff matches the locked pilot decision", {
  expect_equal(AGG_WALD_CUTOFF_DEFAULT, 1)
  expect_equal(AGG_WALD_LAMBDA, 1 / AGG_WALD_CUTOFF_DEFAULT)
  expect_equal(AGG_WALD_CUTOFF, 1 / AGG_WALD_LAMBDA)
})

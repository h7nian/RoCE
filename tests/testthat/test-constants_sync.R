# test-constants_sync.R - Verify R and C++ constants remain in sync
#
# Reads src/numerical_constants.hpp and compares parsed values against
# the corresponding R constants in R/constants.R.
# Prevents silent drift between the two files.

test_that("R and C++ numerical constants are in sync", {
  # Read and parse C++ constants from header file
  hpp_path <- file.path("..", "..", "src", "numerical_constants.hpp")
  if (!file.exists(hpp_path)) {
    hpp_path <- file.path(system.file(package = "FACEC"), "..", "src", "numerical_constants.hpp")
  }
  # Also try the common testthat working directory
  if (!file.exists(hpp_path)) {
    hpp_path <- "src/numerical_constants.hpp"
  }
  skip_if_not(file.exists(hpp_path), "Cannot locate numerical_constants.hpp")
  
  hpp_lines <- readLines(hpp_path)
  
  # Helper: extract constexpr values from C++ header
  parse_cpp_constant <- function(name) {
    pattern <- sprintf("constexpr\\s+\\w+\\s+%s\\s*=\\s*([^;]+);", name)
    match <- grep(pattern, hpp_lines, value = TRUE)
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

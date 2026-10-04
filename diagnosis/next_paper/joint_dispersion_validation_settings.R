make_joint_dispersion_validation_settings <- function() {
  settings <- expand.grid(num_sources = c(4L, 16L, 64L), shared_tate_variance = c(.004, .04),
    profile = c("treated_half", "opposite_halves", "three_quarters", "majority_guarantee",
                "treated_only", "overstated_validity"), local_bias = c(0, .5, 2, 4), stringsAsFactors = FALSE)
  # These configurations duplicate an existing procedure AND its Gaussian data.
  duplicate <- (settings$profile %in% c("treated_only", "overstated_validity") & settings$local_bias == 0) |
    (settings$profile == "majority_guarantee" & settings$num_sources == 4L)
  settings <- settings[!duplicate, ]
  count <- settings$num_sources
  quarter_profiles <- settings$profile %in% c("three_quarters", "majority_guarantee", "overstated_validity")
  settings$valid_mu1_minimum <- ifelse(quarter_profiles, 3L*count/4L, count/2L)
  settings$valid_mu0_minimum <- ifelse(settings$profile == "treated_half", count, settings$valid_mu1_minimum)
  majority <- settings$profile == "majority_guarantee"
  settings$valid_mu1_minimum[majority] <- settings$valid_mu0_minimum[majority] <- floor(count[majority]/2L)+1L
  settings$unshifted_mu1_count <- ifelse(quarter_profiles,3L*count/4L,count/2L)
  settings$unshifted_mu0_count <- ifelse(settings$profile %in% c("treated_half","treated_only"),count,settings$unshifted_mu1_count)
  violation <- settings$profile == "overstated_validity"
  settings$unshifted_mu1_count[violation] <- settings$unshifted_mu0_count[violation] <- count[violation]/4L
  settings$assumption_valid <- settings$local_bias == 0 |
    (settings$unshifted_mu1_count >= settings$valid_mu1_minimum & settings$unshifted_mu0_count >= settings$valid_mu0_minimum)
  rownames(settings) <- NULL
  settings
}

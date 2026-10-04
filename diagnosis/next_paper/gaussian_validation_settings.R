make_gaussian_validation_settings <- function(profile = c("original", "validity_fraction")) {
  profile <- match.arg(profile)
  if (profile == "original") {
    settings <- expand.grid(num_sources = c(4L, 8L, 16L, 32L, 64L),
      shared_correlation = c(0, .25, .5, .9), bias_scale = c(0, .5, 1, 2, 4), stringsAsFactors = FALSE)
    settings$bias_pattern <- "positive"
    alternating <- settings[settings$bias_scale %in% c(1, 2), ]
    alternating$bias_pattern <- "alternating"
    settings <- rbind(settings, alternating)
    settings$unshifted_fraction <- .75
    settings$assumption_valid <- TRUE
    violations <- settings[settings$shared_correlation == .5 & settings$bias_scale == 2 &
                            settings$bias_pattern == "positive", ]
    violations$unshifted_fraction <- .25
    violations$assumption_valid <- FALSE
    settings <- rbind(settings, violations)
    settings$valid_minimum <- ceiling(.75 * settings$num_sources)
    settings$votes_required <- floor(settings$num_sources / 2) + 1L
    settings$seed <- 70100L + seq_len(nrow(settings))
  } else {
    settings <- expand.grid(num_sources = c(4L, 8L, 16L, 32L, 64L),
      shared_correlation = c(0, .5), bias_scale = c(0, .5, 2, 4),
      unshifted_fraction = c(.25, .5, .75), vote_rule = c("partial", "all_valid"),
      stringsAsFactors = FALSE)
    settings$bias_pattern <- "positive"
    settings$assumption_valid <- TRUE
    settings$valid_minimum <- ceiling(settings$unshifted_fraction * settings$num_sources)
    settings$votes_required <- ifelse(settings$vote_rule == "all_valid", settings$valid_minimum,
                                      pmax(1, floor(.6 * settings$valid_minimum)))
    # q=1 gives the same procedure for both vote rules; do not repeat it.
    identity_columns <- setdiff(names(settings), "vote_rule")
    settings <- settings[!duplicated(settings[identity_columns]), ]
    # Different quorum rules evaluate the SAME simulated summaries. At zero
    # bias, different declared valid fractions also share the same data seed.
    data_key <- paste(settings$num_sources, settings$shared_correlation, settings$bias_scale,
      ifelse(settings$bias_scale == 0, settings$num_sources, settings$valid_minimum), sep = ":")
    settings$seed <- 81100L + match(data_key, unique(data_key))
  }
  settings$cell_id <- seq_len(nrow(settings))
  settings$n_per_site <- 1000L
  rownames(settings) <- NULL
  settings
}

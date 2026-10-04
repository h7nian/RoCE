#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if (length(arguments)!=1L) stop("Usage: settings_file")
source(arguments[[1L]])
library(testthat)
test_that("joint settings separate declared validity from actual bias patterns", {
  settings <- make_joint_dispersion_validation_settings()
  expect_equal(nrow(settings),124L)
  expect_setequal(settings$num_sources,c(4L,16L,64L))
  expect_true(all(settings$assumption_valid == (settings$profile != "overstated_validity")))
  majority <- settings[settings$profile=="majority_guarantee",]
  expect_true(all(majority$valid_mu1_minimum==floor(majority$num_sources/2)+1))
  expect_true(all(majority$unshifted_mu1_count==3*majority$num_sources/4))
  single <- settings[settings$profile=="treated_only",]
  expect_true(all(single$unshifted_mu0_count==single$num_sources))
  expect_true(all(single$valid_mu0_minimum==single$num_sources/2))
  identity <- settings[c("num_sources","shared_tate_variance","local_bias",
                         "valid_mu1_minimum","valid_mu0_minimum")]
  identity$shifted1 <- ifelse(settings$local_bias==0,0,settings$num_sources-settings$unshifted_mu1_count)
  identity$shifted0 <- ifelse(settings$local_bias==0,0,settings$num_sources-settings$unshifted_mu0_count)
  expect_false(anyDuplicated(identity)>0)
})
cat("JOINT_DISPERSION_SETTINGS_PASSED\n")

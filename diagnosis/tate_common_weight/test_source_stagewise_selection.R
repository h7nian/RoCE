library(testthat)
root <- if (file.exists("source_stagewise_selection.R")) "." else "diagnosis/tate_common_weight"
source(file.path(root,"target_projection_validation.R"))
source(file.path(root,"source_stagewise_selection.R"))
stagewise_fixture <- function() {
  x <- expand.grid(final_candidate=c("null","c0.5","c1","c2"),
    initial_candidate=c("null","c0.5","c1","c2"),arm=0:1,fold=2:5,stringsAsFactors=FALSE)
  x$n_validation <- c(400,400,398,402)[x$fold-1L]
  x$final_risk <- x$initial_risk <- 0
  x$final_failed <- x$initial_failed <- FALSE
  x
}
test_that("selection follows the chosen downstream branch", {
  x<-stagewise_fixture()
  x$final_risk[x$final_candidate=="c1"] <- -.1
  x$initial_risk[x$final_candidate=="c1" & x$initial_candidate=="c0.5"] <- -.2
  x$initial_risk[x$final_candidate=="c2" & x$initial_candidate=="c2"] <- -100
  s<-.select_source_stagewise_candidates(x)
  expect_identical(s$final_selection$selected,"c1")
  expect_identical(s$initial_selection$selected,"c0.5")
  expect_equal(.select_source_stagewise_candidates(x[nrow(x):1,]),s)
})
test_that("failed upstream candidates do not invalidate final fits", {
  x<-stagewise_fixture(); x$final_risk[x$final_candidate=="c1"] <- -.1
  x$initial_risk[x$final_candidate=="c1" & x$initial_candidate=="c0.5"] <- -1
  x$initial_failed[x$final_candidate=="c1" & x$initial_candidate=="c0.5" & x$fold==2] <- TRUE
  s<-.select_source_stagewise_candidates(x)
  expect_identical(s$final_selection$selected,"c1")
  expect_identical(s$initial_selection$selected,"null")
  expect_false(s$initial_selection$scores$eligible[2])
})
test_that("incomplete and inconsistent candidate tables fail closed", {
  x<-stagewise_fixture()
  expect_identical(.select_source_stagewise_candidates(x)$final_selection$selected,"null")
  expect_error(.select_source_stagewise_candidates(x[-1,]),"complete")
  x$final_risk[2] <- -.1
  expect_error(.select_source_stagewise_candidates(x),"across initial")
  x<-stagewise_fixture(); x$initial_risk[1] <- 1
  expect_error(.select_source_stagewise_candidates(x),"null initial")
})

test_that("one failed final arm invalidates that candidate across folds", {
  x<-stagewise_fixture(); x$final_risk[x$final_candidate=="c1"] <- -1
  x$final_failed[x$final_candidate=="c1" & x$arm==1 & x$fold==3] <- TRUE
  s<-.select_source_stagewise_candidates(x)
  expect_identical(s$final_selection$selected,"null")
  expect_false(s$final_selection$scores$eligible[s$final_selection$scores$candidate=="c1"])
  expect_identical(s$initial_selection$selected,"null")
})

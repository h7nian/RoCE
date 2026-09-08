library(testthat)

.audit_test_repo <- if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else {
  normalizePath(file.path("..", ".."), mustWork = TRUE)
}
setwd(.audit_test_repo)

source("diagnosis/tate_common_weight/audit_full_refit_intermediates.R")
source("scripts/slurm/result_provenance.R")

.audit_test_draw <- file.path(
  "results", "direct_tate_mc500_b5000", "full_refit_calibration_v18",
  "seed_000001_rho_0_draw_0000"
)
.audit_test_positive_draw <- file.path(
  "results", "direct_tate_mc500_b5000", "full_refit_calibration_v18",
  "seed_000001_rho_0_draw_0001"
)
.audit_test_library <- file.path(
  "results", "direct_tate_mc500_b5000",
  "Rlib_weight_bootstrap_20260904_v18"
)

.copy_draw_fixture <- function(source = .audit_test_draw) {
  destination <- tempfile("full_refit_audit_fixture_")
  dir.create(destination)
  files <- list.files(source, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  stopifnot(all(file.copy(files, destination, recursive = TRUE)))
  destination
}

.rehash_draw_fixture <- function(directory) {
  payloads <- c("result.csv", "draw.rds")
  hashes <- vapply(
    file.path(directory, payloads), roce_sha256_file, character(1L)
  )
  writeLines(
    paste(hashes, payloads, sep = "  "),
    file.path(directory, "sha256.txt")
  )
  invisible(directory)
}

.mutate_draw_rds <- function(directory, mutate) {
  path <- file.path(directory, "draw.rds")
  saved <- readRDS(path)
  saved <- mutate(saved)
  saveRDS(saved, path)
  .rehash_draw_fixture(directory)
}

.run_audit_fixture <- function(draw, output = tempfile("full_refit_audit_out_")) {
  audit_full_refit_intermediates_main(c(draw, .audit_test_library, output))
  output
}

test_that("standalone auditor accepts an untouched, rehashed identity fixture", {
  fixture <- .copy_draw_fixture()
  .rehash_draw_fixture(fixture)
  output <- .run_audit_fixture(fixture)

  expect_true(file.exists(file.path(output, "audit_passed.txt")))
  expect_match(
    readLines(file.path(output, "audit_passed.txt"), n = 1L),
    "audit=passed"
  )
})

test_that("standalone auditor rejects validly rehashed wrong row maps and counts", {
  wrong_map <- .copy_draw_fixture()
  .mutate_draw_rds(wrong_map, function(saved) {
    saved$result$original_row_ids$t[[1L]] <- 2L
    saved
  })
  expect_error(.run_audit_fixture(wrong_map), "row maps or multiplicities")

  wrong_count <- .copy_draw_fixture()
  .mutate_draw_rds(wrong_count, function(saved) {
    saved$result$multiplicities$s1[[1L]] <- 2L
    saved
  })
  expect_error(.run_audit_fixture(wrong_count), "row maps or multiplicities")
})

test_that("cross-fold origins fail even when replay and saved maps agree", {
  saved <- readRDS(file.path(.audit_test_positive_draw, "draw.rds"))
  source_saved <- readRDS(file.path(saved$source_bundle, "artifacts.rds"))
  reference <- source_saved$group_result$artifacts[["0"]]$
    direct_tate_results$one_round_crossfit
  info <- reference$intermediates$fold_info
  sources <- names(reference$weights)
  folds <- c(
    list(t = lapply(info, `[[`, "target_idx")),
    setNames(lapply(seq_along(sources), function(j) {
      lapply(info, function(x) x$source_idx[[j]])
    }), sources)
  )
  maps <- saved$result$original_row_ids
  positions_fold_1 <- folds$t[[1L]]
  origin_fold_2 <- folds$t[[2L]][[1L]]
  maps$t[positions_fold_1[[1L]]] <- origin_fold_2
  multiplicities <- lapply(maps, function(x) tabulate(x, nbins = length(x)))
  replay <- list(original_row_ids = maps, multiplicities = multiplicities)

  expect_error(
    roce_audit_full_refit_resampling(maps, multiplicities, replay, folds),
    "crossed its original outer fold"
  )
})

test_that("CSV/RDS scalar disagreements fail after valid checksum rewrite", {
  fixture <- .copy_draw_fixture()
  csv_path <- file.path(fixture, "result.csv")
  csv <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  csv$estimate_reference <- csv$estimate_reference + 0.01
  write.csv(csv, csv_path, row.names = FALSE)
  .rehash_draw_fixture(fixture)

  expect_error(.run_audit_fixture(fixture), "result.csv and draw.rds")
})

test_that("altered raw moments fail independently of valid bundle hashes", {
  fixture <- .copy_draw_fixture()
  .mutate_draw_rds(fixture, function(saved) {
    saved$result$fitted$intermediates$fold_info[[1L]]$V_ot <-
      saved$result$fitted$intermediates$fold_info[[1L]]$V_ot + 0.01
    saved
  })

  expect_error(.run_audit_fixture(fixture), "raw influence moments")
})

test_that("missing, invalid-dimensional, and misaligned source IDs fail closed", {
  missing_id <- .copy_draw_fixture()
  .mutate_draw_rds(missing_id, function(saved) {
    saved$result$fitted$intermediates$fold_info[[1L]]$source_idx[[1L]] <- NULL
    saved
  })
  expect_error(.run_audit_fixture(missing_id), "source_idx|fold structure")

  invalid_dimension <- .copy_draw_fixture()
  .mutate_draw_rds(invalid_dimension, function(saved) {
    saved$result$fitted$intermediates$inner_fold_info[[1L]]$
      components[[1L]]$C_cross <- matrix(0, 1L, 1L)
    saved
  })
  expect_error(
    .run_audit_fixture(invalid_dimension),
    "covariance|dimensions|fold structure"
  )

  misaligned_id <- .copy_draw_fixture()
  .mutate_draw_rds(misaligned_id, function(saved) {
    idx <- saved$result$fitted$intermediates$inner_fold_info[[1L]]$
      components[[1L]]$source_idx[[1L]]
    saved$result$fitted$intermediates$inner_fold_info[[1L]]$
      components[[1L]]$source_idx[[1L]] <- rev(idx)
    saved
  })
  expect_error(.run_audit_fixture(misaligned_id), "exactly match|fold structure")
})

test_that("auditor never overwrites an existing output destination", {
  output <- tempfile("existing_full_refit_audit_")
  dir.create(output)
  marker <- file.path(output, "keep.txt")
  writeLines("preserve", marker)

  expect_error(.run_audit_fixture(.audit_test_draw, output), "already exists")
  expect_identical(readLines(marker), "preserve")
})

test_that("declared identity scope rejects a wrong fitted fold weight", {
  fixture <- .copy_draw_fixture()
  .mutate_draw_rds(fixture, function(saved) {
    saved$result$fitted$fold_weights[1L, 1L] <-
      saved$result$fitted$fold_weights[1L, 1L] + 0.01
    saved
  })

  expect_error(.run_audit_fixture(fixture), "identity draw")
})

test_that("origin-grouped partition records are independently validated", {
  valid_record <- list(
    caller = "synthetic", n_folds = 3L,
    cv_group_id = c(1L, 1L, 2L, 2L, 3L, 3L),
    cv_fold_id = c(1L, 1L, 2L, 2L, 3L, 3L), valid = TRUE
  )
  result <- list(
    nuisance_cv_grouping = "origin",
    nuisance_cv_partitions = list(valid_record)
  )
  checked <- roce_audit_nuisance_cv_partitions(result, "origin_grouped")
  expect_identical(checked$n_records, 1L)
  expect_identical(checked$n_explicit, 1L)

  crossing <- result
  crossing$nuisance_cv_partitions[[1L]]$cv_fold_id[[2L]] <- 2L
  expect_error(
    roce_audit_nuisance_cv_partitions(crossing, "origin_grouped"),
    "independent group/fold validation"
  )

  false_flag <- result
  false_flag$nuisance_cv_partitions[[1L]]$valid <- FALSE
  expect_error(
    roce_audit_nuisance_cv_partitions(false_flag, "origin_grouped"),
    "invalid nuisance CV partition record"
  )

  empty <- result
  empty$nuisance_cv_partitions <- list()
  expect_error(
    roce_audit_nuisance_cv_partitions(empty, "origin_grouped"),
    "retain nonempty"
  )

  missing_explicit <- result
  missing_explicit$nuisance_cv_partitions[[1L]]["cv_fold_id"] <- list(NULL)
  expect_error(
    roce_audit_nuisance_cv_partitions(missing_explicit, "origin_grouped"),
    "independent group/fold validation"
  )

  unique_null <- result
  unique_null$nuisance_cv_partitions[[1L]]$cv_group_id <- 1:6
  unique_null$nuisance_cv_partitions[[1L]]["cv_fold_id"] <- list(NULL)
  checked_unique <- roce_audit_nuisance_cv_partitions(
    unique_null, "origin_grouped"
  )
  expect_identical(checked_unique$n_explicit, 0L)
})

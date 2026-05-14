# test-checkpoint.R - Round-trip and atomicity tests for checkpoint helpers.
#
# Covers (exported functions in R/checkpoint.R):
#   - init_checkpoint_config
#   - save_checkpoint
#   - load_checkpoint
#   - check_preempt_signal
#   - signal_checkpoint_saved
#   - cleanup_checkpoint
#
# These helpers are the load-bearing piece of the SLURM preempt-and-requeue
# flow. A silent corruption (e.g. a partial write or a non-atomic rename) would
# go unnoticed in production until a job actually got preempted, at which point
# the simulation state would be lost. This file pins the contract.

library(testthat)

# Package loaded by helper-load.R.

# ============================================================================
# init_checkpoint_config
# ============================================================================

test_that("init_checkpoint_config creates the directory and builds the expected paths", {
  base_dir <- tempfile(pattern = "facec_ckpt_init_")

  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "demo_setting",
                                   job_id = "job-12345")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  expect_true(dir.exists(base_dir),
              info = "init_checkpoint_config must create the checkpoint directory.")
  expect_equal(config$dir, base_dir)
  expect_equal(config$file,
               file.path(base_dir, "checkpoint_demo_setting_job-12345.RData"))
  expect_equal(config$preempt_signal,
               file.path(base_dir, ".preempt_signal_job-12345"))
  expect_equal(config$saved_signal,
               file.path(base_dir, ".checkpoint_saved_job-12345"))
  expect_true(is.numeric(config$setting_interval) && config$setting_interval > 0)
  expect_true(is.numeric(config$sim_interval) && config$sim_interval > 0)
})

# ============================================================================
# save_checkpoint -> load_checkpoint round-trip
# ============================================================================

test_that("save_checkpoint -> load_checkpoint round-trips a non-trivial state list", {
  base_dir <- tempfile(pattern = "facec_ckpt_roundtrip_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "rt_setting",
                                   job_id = "rt_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  state <- list(
    current_setting_idx = 2L,
    total_settings      = 5L,
    current_sim_idx     = 7L,
    n_sims              = 50L,
    sim_results         = data.frame(method = c("one_round_crossfit",
                                                "two_round_crossfit"),
                                     estimate = c(0.31, 0.29),
                                     stringsAsFactors = FALSE),
    payload             = list(seed = 42L, message = "hello"),
    timestamp           = Sys.time()
  )

  expect_true(save_checkpoint(state, config$file, sim_level = TRUE))
  expect_true(file.exists(config$file))
  expect_false(file.exists(paste0(config$file, ".tmp")),
               info = "save_checkpoint must remove its atomic-write temp file.")

  restored <- load_checkpoint(config$file)
  expect_equal(restored$current_setting_idx, state$current_setting_idx)
  expect_equal(restored$total_settings,      state$total_settings)
  expect_equal(restored$current_sim_idx,     state$current_sim_idx)
  expect_equal(restored$n_sims,              state$n_sims)
  expect_equal(restored$sim_results,         state$sim_results)
  expect_equal(restored$payload,             state$payload)
})

test_that("load_checkpoint returns NULL when no checkpoint file exists", {
  base_dir <- tempfile(pattern = "facec_ckpt_missing_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "missing_setting",
                                   job_id = "missing_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  expect_null(load_checkpoint(config$file))
})

test_that("load_checkpoint fails fast when an existing checkpoint file is corrupt", {
  base_dir <- tempfile(pattern = "facec_ckpt_corrupt_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "corrupt_setting",
                                   job_id = "corrupt_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  # Write a binary blob that is NOT a serialized `state` list. base::load() may
  # emit its own "no readable RDS" warning before raising the error; we
  # silence the warning here so the test focuses on the user-facing failure.
  writeBin(as.raw(c(1, 2, 3, 4)), config$file)
  suppressWarnings(
    expect_error(load_checkpoint(config$file))
  )
})

# ============================================================================
# preempt-signal helpers
# ============================================================================

test_that("check_preempt_signal flips with the on-disk presence of the signal file", {
  base_dir <- tempfile(pattern = "facec_ckpt_preempt_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "preempt_setting",
                                   job_id = "preempt_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  expect_false(check_preempt_signal(config$preempt_signal))
  file.create(config$preempt_signal)
  expect_true(check_preempt_signal(config$preempt_signal))
})

test_that("signal_checkpoint_saved creates the saved-signal file", {
  base_dir <- tempfile(pattern = "facec_ckpt_saved_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "saved_setting",
                                   job_id = "saved_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  expect_false(file.exists(config$saved_signal))
  expect_true(signal_checkpoint_saved(config$saved_signal))
  expect_true(file.exists(config$saved_signal))
})

# ============================================================================
# cleanup_checkpoint
# ============================================================================

test_that("cleanup_checkpoint removes every artifact it knows about", {
  base_dir <- tempfile(pattern = "facec_ckpt_cleanup_")
  config <- init_checkpoint_config(checkpoint_dir = base_dir,
                                   setting_id = "cleanup_setting",
                                   job_id = "cleanup_job")
  on.exit(unlink(base_dir, recursive = TRUE), add = TRUE)

  # Plant all three artifacts.
  save_checkpoint(list(current_setting_idx = 1L, total_settings = 1L,
                       current_sim_idx = NULL, n_sims = NULL,
                       sim_results = NULL),
                  config$file, sim_level = FALSE)
  file.create(config$preempt_signal)
  signal_checkpoint_saved(config$saved_signal)

  expect_true(all(file.exists(c(config$file, config$preempt_signal, config$saved_signal))))

  cleanup_checkpoint(config$file, config$preempt_signal, config$saved_signal)

  expect_false(any(file.exists(c(config$file, config$preempt_signal, config$saved_signal))),
               info = "cleanup_checkpoint must remove every artifact it was given.")
})

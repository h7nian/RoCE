# checkpoint.R - Checkpoint/restart utilities for SLURM preemption support
#
# Provides atomic checkpoint save/load, preemption signal detection, and cleanup.
# All functions take explicit path arguments — no global state dependency.

#' Initialize checkpoint directory and file paths
#'
#' @param checkpoint_dir Base directory for checkpoint files.
#' @param setting_id String identifying the current simulation setting.
#' @param job_id SLURM job ID or unique identifier.
#' @return Named list of checkpoint configuration (dir, file, preempt, saved paths, intervals).
#' @export
init_checkpoint_config <- function(checkpoint_dir = "checkpoints",
                                   setting_id = "", job_id = "") {
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (!dir.exists(checkpoint_dir)) {
    stop(sprintf("init_checkpoint_config: failed to create checkpoint directory '%s'.",
                 checkpoint_dir), call. = FALSE)
  }

  config <- list(
    dir              = checkpoint_dir,
    file             = file.path(checkpoint_dir, sprintf("checkpoint_%s_%s.RData", setting_id, job_id)),
    preempt_signal   = file.path(checkpoint_dir, sprintf(".preempt_signal_%s", job_id)),
    saved_signal     = file.path(checkpoint_dir, sprintf(".checkpoint_saved_%s", job_id)),
    setting_interval = CHECKPOINT_SETTING_INTERVAL,   # Save checkpoint every N parameter settings
    sim_interval     = CHECKPOINT_SIM_INTERVAL   # Save checkpoint every N simulations within a setting
  )
  return(config)
}

#' Save checkpoint state (atomic: write to temp then rename)
#'
#' @param state List containing current simulation state.
#' @param checkpoint_file Path to checkpoint RData file.
#' @param sim_level Logical; if TRUE, print a compact simulation-level message.
#' @return TRUE invisibly on success; throws an error on failure.
#' @export
save_checkpoint <- function(state, checkpoint_file, sim_level = FALSE) {
  temp_file <- paste0(checkpoint_file, ".tmp")
  save(state, file = temp_file)
  renamed <- file.rename(temp_file, checkpoint_file)
  if (!isTRUE(renamed)) {
    stop(sprintf("save_checkpoint: failed to move temporary checkpoint '%s' to '%s'.",
                 temp_file, checkpoint_file), call. = FALSE)
  }

  if (sim_level) {
    cat(sprintf("  Checkpoint saved: setting %d, simulation %d/%d\n",
                state$current_setting_idx, state$current_sim_idx, state$n_sims))
  } else {
    cat(sprintf("Checkpoint saved: %s (setting %d/%d completed)\n",
                checkpoint_file, state$current_setting_idx, state$total_settings))
  }
  invisible(TRUE)
}

#' Load checkpoint state if the file exists
#'
#' @param checkpoint_file Path to checkpoint RData file.
#' @return Restored state list, or NULL if no checkpoint found. Throws an
#'   error if an existing checkpoint cannot be loaded.
#' @export
load_checkpoint <- function(checkpoint_file) {
  if (!file.exists(checkpoint_file)) return(NULL)

  checkpoint_env <- new.env(parent = emptyenv())
  tryCatch(
    load(checkpoint_file, envir = checkpoint_env),  # restores 'state'
    warning = function(w) {
      stop(sprintf("load_checkpoint: failed to load checkpoint '%s': %s",
                   checkpoint_file, conditionMessage(w)), call. = FALSE)
    },
    error = function(e) {
      stop(sprintf("load_checkpoint: failed to load checkpoint '%s': %s",
                   checkpoint_file, conditionMessage(e)), call. = FALSE)
    }
  )
  if (!exists("state", envir = checkpoint_env, inherits = FALSE) ||
      !is.list(checkpoint_env$state)) {
    stop(sprintf("load_checkpoint: '%s' did not contain a valid checkpoint state.",
                 checkpoint_file), call. = FALSE)
  }

  state <- checkpoint_env$state
  if (!is.null(state$sim_results) && !is.null(state$current_sim_idx)) {
    cat(sprintf("Checkpoint loaded: %s\n", checkpoint_file))
    cat(sprintf("   Resuming setting %d/%d from simulation %d/%d\n",
                state$current_setting_idx, state$total_settings,
                state$current_sim_idx + 1, state$n_sims))
  } else {
    cat(sprintf("Checkpoint loaded: %s (resuming from setting %d/%d)\n",
                checkpoint_file, state$current_setting_idx + 1, state$total_settings))
  }
  return(state)
}

#' Check if a SLURM preemption signal was received
#'
#' @param preempt_signal_file Path to the preemption signal file.
#' @return TRUE if the signal file exists.
#' @export
check_preempt_signal <- function(preempt_signal_file) {
  return(file.exists(preempt_signal_file))
}

#' Signal that a checkpoint has been saved (for shell-script coordination)
#'
#' @param saved_signal_file Path to the saved-signal file.
#' @export
signal_checkpoint_saved <- function(saved_signal_file) {
  created <- file.create(saved_signal_file)
  if (!isTRUE(created)) {
    stop(sprintf("signal_checkpoint_saved: failed to create '%s'.",
                 saved_signal_file), call. = FALSE)
  }
  invisible(TRUE)
}

#' Clean up all checkpoint artifacts after successful completion
#'
#' @param checkpoint_file Path to checkpoint RData file.
#' @param preempt_signal_file Path to preemption signal file.
#' @param saved_signal_file Path to saved-signal file.
#' @export
cleanup_checkpoint <- function(checkpoint_file, preempt_signal_file, saved_signal_file) {
  for (f in c(checkpoint_file, preempt_signal_file, saved_signal_file)) {
    if (file.exists(f)) file.remove(f)
  }
  cat("Checkpoint files removed after successful completion.\n")
}

# Shared atomic-output helper for standalone Slurm R drivers.
#
# Keeping the cleanup handler inside this function is intentional: on.exit()
# registered at the top level of an Rscript is not run when the script exits.
roce_write_atomic_directory <- function(output_directory, writer,
                                          caller = "atomic output") {
  if (!is.character(output_directory) || length(output_directory) != 1L ||
      is.na(output_directory) || !nzchar(output_directory)) {
    stop(caller, ": output_directory must be one nonempty path.",
         call. = FALSE)
  }
  if (!is.function(writer)) {
    stop(caller, ": writer must be a function.", call. = FALSE)
  }
  if (dir.exists(output_directory) || file.exists(output_directory)) {
    stop(caller, ": output already exists: ", output_directory,
         call. = FALSE)
  }

  parent_directory <- dirname(output_directory)
  dir.create(parent_directory, recursive = TRUE, showWarnings = FALSE)
  staging_directory <- tempfile(
    pattern = paste0(".", basename(output_directory), "_"),
    tmpdir = parent_directory
  )
  if (!dir.create(staging_directory, showWarnings = FALSE)) {
    stop(caller, ": failed to create staging directory.", call. = FALSE)
  }
  on.exit({
    if (dir.exists(staging_directory)) {
      unlink(staging_directory, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  writer(staging_directory)
  staged_files <- list.files(
    staging_directory, all.files = TRUE, no.. = TRUE
  )
  if (length(staged_files) == 0L) {
    stop(caller, ": writer produced an empty staging directory.",
         call. = FALSE)
  }
  if (!file.rename(staging_directory, output_directory)) {
    stop(caller, ": failed to atomically install output directory: ",
         output_directory, call. = FALSE)
  }
  invisible(output_directory)
}

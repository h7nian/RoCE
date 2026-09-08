#!/usr/bin/env Rscript

# Isolated integration tests for the audit-array scheduling wrapper. Shell
# module/Rscript calls are stubbed; no RoCE model or scientific audit is run.
main <- function() {
  project <- normalizePath(".", mustWork = TRUE)
  wrapper <- file.path(
    project, "diagnosis/tate_common_weight",
    "run_independent_inference_pilot_audit_array.sh"
  )
  if (!file.exists(wrapper)) {
    stop("audit-array wrapper is not available yet: ", wrapper, call. = FALSE)
  }
  fixture <- tempfile("roce_inference_audit_array_")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  write_fixture <- function(path, lines = "fixture") {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
  }
  fixture_project <- file.path(fixture, "project")
  dir.create(file.path(fixture_project, "scripts", "slurm"),
             recursive = TRUE)
  stopifnot(file.copy(
    file.path(project, "scripts/slurm/package_library_utils.sh"),
    file.path(fixture_project, "scripts/slurm/package_library_utils.sh")
  ))
  audit_script <- file.path(
    fixture_project, "diagnosis/tate_common_weight",
    "audit_independent_inference_pilot.R"
  )
  write_fixture(audit_script, "stop('the Rscript stub should intercept this')")

  pilot_root <- file.path(fixture, "pilot")
  dir.create(pilot_root)
  write_fixture(file.path(pilot_root, "manifest.csv"), "task_id")
  for (task in c(1L, 11L, 25L, 100L)) {
    dir.create(file.path(
      pilot_root, sprintf("seed_%06d", 10000L + task)
    ))
  }
  library <- file.path(fixture, "library")
  package_files <- c(
    "DESCRIPTION", "libs/RoCE.so", "R/RoCE.rdb", "R/RoCE.rdx"
  )
  for (path in package_files) {
    write_fixture(file.path(library, "RoCE", path), path)
  }

  call_log <- file.path(fixture, "rscript_calls.txt")
  module_log <- file.path(fixture, "module_calls.txt")
  stubs <- file.path(fixture, "shell_stubs.sh")
  write_fixture(stubs, c(
    "module() { printf '%s\\n' \"$*\" >> \"${STUB_MODULE_LOG}\"; }",
    paste0("Rscript() { printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' ",
      "\"$1\" \"$2\" \"$3\" \"$4\" \"${OMP_NUM_THREADS:-}\" ",
      "\"${OPENBLAS_NUM_THREADS:-}\" \"${MKL_NUM_THREADS:-}\" ",
      ">> \"${STUB_RSCRIPT_LOG}\"; return \"${STUB_RSCRIPT_STATUS:-0}\"; }"
    )
  ))

  case_count <- 0L
  run_case <- function(label, task_id, root = pilot_root,
                       package_library = library, expected_status = 0L,
                       expect_call = expected_status == 0L,
                       prepare = NULL, expected_bundle = NULL,
                       expected_output = NULL) {
    case_count <<- case_count + 1L
    if (file.exists(call_log)) unlink(call_log)
    if (file.exists(module_log)) unlink(module_log)
    if (is.function(prepare)) prepare()
    environment <- c(
      BASH_ENV = stubs,
      ROCE_PROJECT_ROOT = fixture_project,
      SLURM_ARRAY_TASK_ID = task_id,
      STUB_RSCRIPT_LOG = call_log,
      STUB_MODULE_LOG = module_log,
      STUB_RSCRIPT_STATUS = as.character(expected_status)
    )
    output <- suppressWarnings(system2(
      "bash", c(shQuote(wrapper), shQuote(root), shQuote(package_library)),
      env = paste(names(environment), shQuote(environment), sep = "="),
      stdout = TRUE, stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    calls <- if (file.exists(call_log)) readLines(call_log, warn = FALSE) else character()
    if (status != expected_status || (length(calls) == 1L) != expect_call) {
      stop(label, ": unexpected status or Rscript call\n",
           paste(output, collapse = "\n"), call. = FALSE)
    }
    if (expect_call) {
      fields <- strsplit(calls[[1L]], "\t", fixed = TRUE)[[1L]]
      if (length(fields) != 7L ||
          !identical(fields[[1L]],
                     "diagnosis/tate_common_weight/audit_independent_inference_pilot.R") ||
          !identical(normalizePath(fields[[2L]], mustWork = TRUE),
                     normalizePath(expected_bundle, mustWork = TRUE)) ||
          !identical(normalizePath(fields[[3L]], mustWork = TRUE),
                     normalizePath(library, mustWork = TRUE)) ||
          !identical(fields[[4L]], expected_output) ||
          !identical(fields[5:7], rep("1", 3L))) {
        stop(label, ": Rscript arguments were not mapped exactly.",
             call. = FALSE)
      }
      modules <- readLines(module_log, warn = FALSE)
      if (!identical(modules, "load R/4.2.2-gcc-8.2.0-vp7tyde")) {
        stop(label, ": fixed R module was not requested exactly.",
             call. = FALSE)
      }
      if (dir.exists(expected_output) || file.exists(expected_output)) {
        stop(label, ": wrapper/stub unexpectedly published audit output.",
             call. = FALSE)
      }
    }
    message("[passed] ", label)
  }

  for (task in c(1L, 11L, 25L, 100L)) {
    seed <- 10000L + task
    bundle <- file.path(pilot_root, sprintf("seed_%06d", seed))
    destination <- file.path(
      pilot_root, "audits", sprintf("seed_%06d", seed)
    )
    run_case(
      paste("positive task mapping", task), as.character(task),
      expected_bundle = bundle, expected_output = destination
    )
  }

  for (invalid in c("0", "01", "1.0", "-1", "101", "text", "")) {
    run_case(paste0("reject invalid array ID '", invalid, "'"), invalid,
             expected_status = 1L, expect_call = FALSE)
  }

  absent_root <- file.path(fixture, "absent_root")
  run_case("reject missing pilot root", "1", root = absent_root,
           expected_status = 1L, expect_call = FALSE)
  no_manifest <- file.path(fixture, "no_manifest")
  dir.create(no_manifest)
  run_case("reject missing manifest", "1", root = no_manifest,
           expected_status = 1L, expect_call = FALSE)
  run_case("reject missing source bundle", "2", expected_status = 1L,
           expect_call = FALSE)
  run_case("reject missing package library", "1",
           package_library = file.path(fixture, "absent_library"),
           expected_status = 1L, expect_call = FALSE)

  incomplete <- file.path(fixture, "incomplete_library")
  for (path in package_files[-4L]) {
    write_fixture(file.path(incomplete, "RoCE", path), path)
  }
  run_case("reject incomplete package library", "1",
           package_library = incomplete, expected_status = 1L,
           expect_call = FALSE)

  existing_output <- file.path(pilot_root, "audits", "seed_010011")
  run_case("reject existing audit output", "11", expected_status = 1L,
           expect_call = FALSE,
           prepare = function() dir.create(existing_output, recursive = TRUE))
  unlink(existing_output, recursive = TRUE)
  existing_lock <- paste0(existing_output, ".lock")
  run_case("reject existing audit lock", "11", expected_status = 1L,
           expect_call = FALSE,
           prepare = function() write_fixture(existing_lock, "locked"))
  unlink(existing_lock)

  run_case(
    "propagate Rscript nonzero status without publishing success", "25",
    expected_status = 37L, expect_call = TRUE,
    expected_bundle = file.path(pilot_root, "seed_010025"),
    expected_output = file.path(pilot_root, "audits", "seed_010025")
  )
  if (dir.exists(file.path(pilot_root, "audits", "seed_010025"))) {
    stop("nonzero Rscript stub unexpectedly published an audit directory.",
         call. = FALSE)
  }
  message("[passed] ", case_count,
          " isolated audit-array wrapper cases; no scientific audit run")
}

main()

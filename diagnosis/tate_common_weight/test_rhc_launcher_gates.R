#!/usr/bin/env Rscript

# Isolated launcher integration tests. The RHC model call and module loading
# are shell stubs; all gate checks and fingerprint functions are real.
main <- function() {
  project <- normalizePath(".", mustWork = TRUE)
  root <- tempfile("roce_rhc_launcher_gates_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_fixture <- function(path, lines) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path)
  }
  shell_output <- function(command) {
    output <- system2("bash", c("-c", shQuote(command)), stdout = TRUE,
                      stderr = TRUE)
    if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
    output
  }
  source_paths <- c("DESCRIPTION", "LICENSE", "NAMESPACE", "README.md",
                    "R/fixture.R", "src/fixture.cpp", "inst/fixture.txt",
                    "man/fixture.Rd", "tests/testthat/fixture.R")
  for (path in source_paths) write_fixture(file.path(root, path), path)
  scripts <- c("run_rhc_direct_tate.sh", "run_rhc_direct_tate.R",
               "direct_tate_task_helpers.R", "result_provenance.R",
               "resource_topology.R", "package_library_utils.sh",
               "run_package_audit_tests.sh", "run_r_cmd_check.sh")
  dir.create(file.path(root, "scripts", "slurm"), recursive = TRUE)
  stopifnot(all(file.copy(file.path(project, "scripts", "slurm", scripts),
                         file.path(root, "scripts", "slurm", scripts))))
  library <- file.path(root, "library")
  package_files <- c("DESCRIPTION", "libs/RoCE.so", "R/RoCE.rdb",
                     "R/RoCE.rdx", "extdata/rhc.csv")
  for (path in package_files) write_fixture(file.path(library, "RoCE", path), path)
  utilities <- file.path(root, "scripts", "slurm", "package_library_utils.sh")
  package_hash <- shell_output(paste("source", shQuote(utilities),
    "; roce_package_fingerprint", shQuote(library)))
  source_hash <- shell_output(paste("source", shQuote(utilities),
    "; roce_package_source_fingerprint", shQuote(root)))
  test_hash <- shell_output(paste("find", shQuote(file.path(root, "tests")),
    "-type f -name '*.R' -print0 | sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}'"))
  test_gate <- file.path(library, "audit_tests_passed.txt")
  check_gate <- file.path(root, "check", "r_cmd_check_passed.txt")
  test_values <- c(
    package_tests = "passed", package_fingerprint = package_hash,
    package_source_fingerprint = source_hash, test_suite_fingerprint = test_hash,
    test_driver_md5 = unname(tools::md5sum(file.path(root, "scripts", "slurm",
                                                   "run_package_audit_tests.sh"))))
  check_values <- c(
    r_cmd_check = "passed", package_source_fingerprint = source_hash,
    check_driver_md5 = unname(tools::md5sum(file.path(root, "scripts", "slurm",
                                                    "run_r_cmd_check.sh"))))
  write_gate <- function(path, values) {
    write_fixture(path, paste(names(values), values, sep = "="))
  }
  reset_gates <- function() {
    write_gate(test_gate, test_values)
    write_gate(check_gate, check_values)
  }
  stub <- file.path(root, "shell_stubs.sh")
  write_fixture(stub, c("module() { :; }",
    "Rscript() { printf '%s\\n' '[stub] RHC execution reached'; }"))
  launcher <- file.path(root, "scripts", "slurm", "run_rhc_direct_tate.sh")
  case_count <- 0L
  run_case <- function(label, expected_success, check_path = check_gate,
                       expected_error = NULL) {
    case_count <<- case_count + 1L
    output_root <- file.path(root, paste0("output_", case_count))
    environment <- c(BASH_ENV = stub, ROCE_PROJECT_ROOT = root,
      ROCE_PROJECT_LIB = library, ROCE_OUTPUT_ROOT = output_root,
      ROCE_PACKAGE_CHECK_GATE = check_path, SLURM_CPUS_PER_TASK = "30")
    output <- suppressWarnings(system2("bash", shQuote(launcher),
      env = paste(names(environment), shQuote(environment), sep = "="),
      stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    reached_model_stub <- any(output == "[stub] RHC execution reached")
    if (!identical(status == 0L, expected_success) ||
        !identical(reached_model_stub, expected_success) ||
        !identical(dir.exists(output_root), expected_success) ||
        (!is.null(expected_error) &&
         !any(grepl(expected_error, output, fixed = TRUE)))) {
      stop(label, ": unexpected launcher outcome\n", paste(output, collapse = "\n"))
    }
    message("[passed] ", label)
  }

  reset_gates()
  run_case("matching gates reach only the model stub", TRUE)
  for (key in names(test_values)) {
    reset_gates()
    bad_values <- test_values; bad_values[[key]] <- "wrong"
    write_gate(test_gate, bad_values)
    run_case(paste("reject mismatched test gate", key), FALSE)
  }
  for (key in names(check_values)) {
    reset_gates()
    bad_values <- check_values; bad_values[[key]] <- "wrong"
    write_gate(check_gate, bad_values)
    run_case(paste("reject mismatched check gate", key), FALSE)
  }
  for (key in names(test_values)) {
    reset_gates()
    write_gate(test_gate, test_values[names(test_values) != key])
    run_case(paste("reject missing test-gate key", key), FALSE,
             expected_error = paste("field", key))
  }
  for (key in names(check_values)) {
    reset_gates()
    write_gate(check_gate, check_values[names(check_values) != key])
    run_case(paste("reject missing check-gate key", key), FALSE,
             expected_error = paste("field", key))
  }
  reset_gates()
  run_case("reject missing check file", FALSE, file.path(root, "absent_gate"))
  run_case("reject empty required check path", FALSE, "")
  write_fixture(test_gate, c(paste(names(test_values), test_values, sep = "="),
                            "package_tests=passed"))
  run_case("reject duplicate matching gate field", FALSE)
  reset_gates()
  write_fixture(check_gate, c(paste(names(check_values), check_values, sep = "="),
                             "r_cmd_check=passed"))
  run_case("reject duplicate matching check field", FALSE,
           expected_error = "field r_cmd_check")
  reset_gates()
  write_fixture(test_gate, c(paste(names(test_values), test_values, sep = "="),
                            "package_fingerprint=wrong"))
  run_case("reject conflicting duplicate package fingerprint", FALSE,
           expected_error = "field package_fingerprint")
  reset_gates()
  unlink(test_gate)
  run_case("reject missing installed-test gate", FALSE)
  reset_gates()
  write_fixture(file.path(library, "RoCE", "libs", "RoCE.so"), "changed bytes")
  run_case("reject changed installed package bytes", FALSE)
  write_fixture(file.path(library, "RoCE", "libs", "RoCE.so"), "libs/RoCE.so")
  write_fixture(file.path(root, "R", "fixture.R"), "changed source")
  run_case("reject changed source bytes", FALSE)
  write_fixture(file.path(root, "R", "fixture.R"), "R/fixture.R")
  write_fixture(file.path(root, "tests", "testthat", "fixture.R"), "changed test")
  changed_source_hash <- shell_output(paste("source", shQuote(utilities),
    "; roce_package_source_fingerprint", shQuote(root)))
  updated_test_values <- test_values
  updated_test_values[["package_source_fingerprint"]] <- changed_source_hash
  write_gate(test_gate, updated_test_values)
  run_case("reject changed tests even with updated source field", FALSE,
           expected_error = "field test_suite_fingerprint")
  write_fixture(file.path(root, "tests", "testthat", "fixture.R"),
                "tests/testthat/fixture.R")
  reset_gates()
  run_case("restored exact bytes pass", TRUE)
  message("[passed] ", case_count, " isolated RHC launcher cases; no model fitted")
}

main()

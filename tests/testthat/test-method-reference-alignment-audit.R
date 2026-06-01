.method_reference_audit_enabled <- function() {
  Sys.getenv("FACEHD_METHOD_REFERENCE_AUDIT", "0") %in%
    c("1", "TRUE", "true", "True")
}

.extract_pdf_text_for_audit <- function(repo_root, ...) {
  pdf_path <- file.path(repo_root, ...)
  extracted_path <- sub("\\.pdf$", ".extracted.txt", pdf_path, ignore.case = TRUE)
  if (file.exists(extracted_path)) {
    return(paste(readLines(extracted_path, warn = FALSE), collapse = "\n"))
  }
  gs <- Sys.which("gs")
  if (!nzchar(gs) && file.exists("/usr/bin/gs")) gs <- "/usr/bin/gs"
  if (!nzchar(gs)) return("")
  out <- tryCatch(
    system2(
      gs,
      c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
        "-sOutputFile=-", pdf_path),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character()
  )
  paste(out, collapse = "\n")
}

test_that("local FACE/SMMAL references contain the audited CV and cross-fitting ideas", {
  skip_if_not(
    .method_reference_audit_enabled(),
    "Set FACEHD_METHOD_REFERENCE_AUDIT=1 to run the method/reference audit."
  )

  repo_root <- normalizePath(file.path(dirname(test_path()), "..", ".."), mustWork = TRUE)
  nospace <- function(x) gsub("[[:space:]]+", "", x)

  face_pdf <- nospace(.extract_pdf_text_for_audit(repo_root, "docs", "FACE.pdf"))
  smmal_pdf <- nospace(.extract_pdf_text_for_audit(repo_root, "docs", "SMMAL.pdf"))
  skip_if_not(nzchar(face_pdf), message = "FACE.pdf text extraction unavailable on this node.")
  skip_if_not(nzchar(smmal_pdf), message = "SMMAL.pdf text extraction unavailable on this node.")

  expect_match(face_pdf, "trainingandvalidationdatasets|samplesplitting", ignore.case = TRUE)
  expect_match(face_pdf, "summarystatistics.*validationdatasets|validationdatasets.*summarystatistics", ignore.case = TRUE)
  expect_match(smmal_pdf, "two[-]*(level|layer)cross|two(level|layer)cross", ignore.case = TRUE)
  expect_match(smmal_pdf, "out-of-two-fold|outoftwofold", ignore.case = TRUE)
  expect_match(smmal_pdf, "calibratedloss|calibratedlosses", ignore.case = TRUE)
})

test_that("main.tex and implementation encode FACE/SMMAL-style alignment", {
  skip_if_not(
    .method_reference_audit_enabled(),
    "Set FACEHD_METHOD_REFERENCE_AUDIT=1 to run the method/reference audit."
  )

  repo_root <- normalizePath(file.path(dirname(test_path()), "..", ".."), mustWork = TRUE)
  read_repo <- function(...) {
    paste(readLines(file.path(repo_root, ...), warn = FALSE), collapse = "\n")
  }

  main_tex <- read_repo("docs", "main.tex")
  agg_code <- read_repo("R", "cross_fitting_aggregation.R")
  cf_code <- read_repo("R", "cross_fitting_algorithms.R")

  expect_match(main_tex, "FACE-style validation for the aggregation penalty", fixed = TRUE)
  expect_match(main_tex, "Fold-summed SMMAL-style calibration", fixed = TRUE)
  expect_match(main_tex, "T_\\{k_1\\}=\\mathcal\\{D\\}\\\\setminus\\mathcal\\{D\\}_\\{k_1\\}|inner-training|inner validation fold")

  expect_match(agg_code, "select_aggregation_lambda_inner_cv", fixed = TRUE)
  expect_match(agg_code, "components\\[-m\\]")
  expect_match(agg_code, "components\\[\\[m\\]\\]")
  expect_match(agg_code, "n_val \\* var_term")

  expect_match(cf_code, "secondary_folds <- setdiff\\(1:n_folds, k1\\)")
  expect_match(cf_code, "training_folds <- setdiff\\(1:n_folds, c\\(k1, k2\\)\\)")
  expect_match(cf_code, ".make_plugin_block_design", fixed = TRUE)
  expect_match(cf_code, "fit_unified_density_ratio", fixed = TRUE)
  expect_match(cf_code, "fit_unified_outcome", fixed = TRUE)
  expect_match(cf_code, "calibrated = TRUE", fixed = TRUE)
})

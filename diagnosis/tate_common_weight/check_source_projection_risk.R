#!/usr/bin/env Rscript
# Verify conditional block losses without asserting a joint scalar potential.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: check_source_projection_risk.R TREATED_PROBE CONTROL_PROBE OUTPUT")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  reports <- list()
  for (arm in 1:0) {
    root <- args[if (arm == 1L) 1L else 2L]
    lines <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(lines[substring(lines, 67L) == "matrix_probe.rds"], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, "matrix_probe.rds")) == expected)
    probe <- readRDS(file.path(root, "matrix_probe.rds"))
    J <- probe$evaluated$jacobian; D <- probe$evaluated$score_gradient
    indices <- probe$state$indices
    block_sign <- ifelse(grepl("^alpha_", names(indices)), -1, 1)
    signs <- numeric(length(D))
    for (k in seq_along(indices)) signs[indices[[k]]] <- block_sign[k]
    # The signed adjoint equation is B a = S D. B need not be symmetric.
    B <- sweep(t(J), 1L, signs, "*")
    asymmetry <- max(abs(B-t(B)))
    for (attempt in probe$attempts) {
      stopifnot(is.null(attempt$failure))
      a <- attempt$coefficients
      residual <- drop(crossprod(J, a))-D
      conditional_error <- finite_difference_error <- 0
      for (k in seq_along(indices)) {
        index <- indices[[k]]; other <- setdiff(seq_along(D), index)
        sign <- block_sign[k]
        H <- sign*J[index, index, drop = FALSE]
        rhs <- sign*(D[index]-drop(crossprod(J[other, index, drop = FALSE], a[other])))
        gradient <- drop(H %*% a[index])-rhs
        conditional_error <- max(conditional_error, abs(gradient-sign*residual[index]))
        direction <- sin(seq_along(index)); direction <- direction/sqrt(sum(direction^2))
        risk <- function(x) .5*sum(x*drop(H %*% x))-sum(rhs*x)
        epsilon <- 1e-5
        numerical <- (risk(a[index]+epsilon*direction)-risk(a[index]-epsilon*direction))/(2*epsilon)
        finite_difference_error <- max(finite_difference_error, abs(numerical-sum(direction*gradient)))
      }
      # A naive joint quadratic necessarily differentiates the symmetric part.
      naive_gradient <- drop((B+t(B)) %*% a)/2-signs*D
      signed_adjoint <- signs*residual
      naive_error <- max(abs(naive_gradient-signed_adjoint))
      stopifnot(conditional_error < 1e-10, finite_difference_error < 1e-7)
      reports[[length(reports)+1L]] <- data.frame(arm, fraction = attempt$summary$fraction,
        signed_adjoint_asymmetry = asymmetry, conditional_gradient_error = conditional_error,
        conditional_directional_error = finite_difference_error,
        naive_joint_quadratic_gradient_error = naive_error)
    }
  }
  reports <- do.call(rbind, reports)
  roce_write_atomic_directory(args[3], function(stage) {
    write.csv(reports, file.path(stage, "risk_gradient_checks.csv"), row.names = FALSE)
    writeLines(c("source_validation_policy_selected=FALSE", "inference_validated=FALSE",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/check_source_projection_risk.R")),
      paste0("input_", 1:2, "_manifest_sha256=", vapply(args[1:2], function(root)
        roce_sha256_file(file.path(root, "sha256.txt")), ""))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "source projection conditional risk check")
  print(reports, row.names = FALSE)
}
if (sys.nframe() == 0L) main()

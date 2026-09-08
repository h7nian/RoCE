#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: check_source_validation_state.R LIBRARY SEED AUDITED_EQUATIONS")
  .libPaths(c(args[1], .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  source("diagnosis/tate_common_weight/source_projection_validation_state.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", equations)
  holdout <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R", holdout)
  artifact <- readRDS(file.path(args[2], "artifacts.rds"))$group_result$artifacts[["0"]]
  probes <- readRDS(file.path(args[3], "validation_equation_states.rds"))
  info <- artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info[[2L]]
  for (arm in 0:1) {
    probe <- probes[[paste0("mu", arm)]]; state <- probe$state
    stopifnot(identical(equations$.evaluate_saved_projection(state, state$theta, TRUE), probe$evaluated))
    target <- equations$.saved_projection_block(artifact$data_split$t, info$target_idx, arm)
    source <- equations$.saved_projection_block(artifact$data_split$s1, info$source_idx[[1L]], arm)
    validation <- .source_projection_validation_state(state, target, source)
    value <- equations$.evaluate_saved_projection(validation, state$theta, TRUE)
    errors <- vapply(1:6, function(k) {
      a <- sin(seq_along(state$theta)*k)/100
      target_score <- holdout$.projection_holdout_scores(state, target, "target", a)
      source_score <- holdout$.projection_holdout_scores(state, source, "source", a)
      moment_error <- abs(sum(value$moments*a)-sum(colMeans(target_score$terms))-sum(colMeans(source_score$terms)))
      direction <- a/sqrt(sum(a^2)); epsilon <- 1e-6
      plus <- equations$.evaluate_saved_projection(validation, state$theta+epsilon*direction)
      minus <- equations$.evaluate_saved_projection(validation, state$theta-epsilon*direction)
      c(moment_error = moment_error,
        jacobian_error = max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(value$jacobian %*% direction))),
        score_gradient_error = abs((plus$score-minus$score)/(2*epsilon)-sum(value$score_gradient*direction)))
    }, numeric(3L))
    stopifnot(max(errors) < 1e-7)
    # Reusing training observations must be rejected before scoring.
    stopifnot(inherits(try(.source_projection_validation_state(state, state$target, source), silent = TRUE), "try-error"))
    cat("arm", arm, "default identity PASS; held-out equation maximum errors:", apply(errors,1,max), "\n")
  }
}
if (sys.nframe() == 0L) main()

#!/usr/bin/env Rscript
# First-split candidate evaluation; intentionally does not select a scale.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(4L,5L)) stop("usage: evaluate_source_projection_candidates.R LIBRARY SEED AUDITED_EQUATIONS OUTPUT [CONFIG]")
  working_config <- if(length(args)==5L) args[5] else "C1"
  .libPaths(c(args[1], .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  for (file in c("sparse_moment_projection.R", "source_final_projection.R", "source_joint_projection.R",
                 "source_projection_validation_state.R")) source(file.path("diagnosis/tate_common_weight", file))
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", equations)
  holdout <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R", holdout)
  checked <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  artifact <- checked(args[2], "artifacts.rds")$group_result$artifacts[["0"]]
  artifact$data_split <- .source_working_basis(artifact$data_split,working_config)
  probes <- checked(args[3], "validation_equation_states.rds")
  info <- artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info[[2L]]
  records <- reports <- risks <- list()
  for (arm in 0:1) {
    probe <- probes[[paste0("mu", arm)]]; state <- probe$state
    target <- equations$.saved_projection_block(artifact$data_split$t, info$target_idx, arm)
    source <- equations$.saved_projection_block(artifact$data_split$s1, info$source_idx[[1L]], arm)
    validation <- .source_projection_validation_state(state, target, source)
    value <- equations$.evaluate_saved_projection(validation, state$theta, TRUE)
    for (scale in c(.5, 1, 2)) {
      candidate <- tryCatch(.fit_source_joint_projection(state, probe$evaluated, scale),
        error = function(e) { if (!inherits(e, "roce_projection_error")) stop(e); e })
      failed <- inherits(candidate, "error")
      key <- paste(arm, scale, sep = ":")
      if (failed) {
        records[[key]] <- list(failure = candidate)
        reports[[key]] <- data.frame(arm, penalty_scale = scale, failed = TRUE,
          original_mean = NA_real_, corrected_mean = NA_real_, target_score_sd = NA_real_, source_score_sd = NA_real_)
        next
      }
      a <- candidate$coefficients
      target_scores <- holdout$.projection_holdout_scores(state, target, "target", a)
      source_scores <- holdout$.projection_holdout_scores(state, source, "source", a)
      corrected <- mean(target_scores$corrected)+mean(source_scores$corrected)
      stopifnot(abs(corrected-value$score+sum(value$moments*a)) < 1e-10)
      for (name in names(state$indices)) {
        index <- state$indices[[name]]; other <- setdiff(seq_along(a), index)
        sign <- if (grepl("^alpha_", name)) -1 else 1
        H <- sign*value$jacobian[index, index, drop = FALSE]
        rhs <- sign*(value$score_gradient[index]-drop(crossprod(value$jacobian[other,index,drop=FALSE],a[other])))
        risk <- .5*sum(a[index]*drop(H %*% a[index]))-sum(rhs*a[index])
        risks[[length(risks)+1L]] <- data.frame(arm, penalty_scale = scale, block = name,
          conditional_risk = risk, downstream_coefficients_vary_across_scales = !name %in% c("alpha_final","gamma_final"))
      }
      records[[key]] <- list(candidate = candidate, target = target_scores, source = source_scores,
        target_ids = target$index, source_ids = source$index)
      reports[[key]] <- data.frame(arm, penalty_scale = scale, failed = FALSE,
        original_mean = value$score, corrected_mean = corrected,
        target_score_sd = sd(target_scores$corrected), source_score_sd = sd(source_scores$corrected))
    }
  }
  reports <- do.call(rbind, reports)
  roce_write_atomic_directory(args[4], function(stage) {
    saveRDS(records, file.path(stage, "candidate_states.rds"))
    write.csv(reports, file.path(stage, "candidate_scores.csv"), row.names = FALSE)
    write.csv(do.call(rbind, risks), file.path(stage, "conditional_block_risks.csv"), row.names = FALSE)
    writeLines(c("penalty_policy_selected=FALSE", "inference_validated=FALSE", "validation_fold=2",paste0("working_config=",working_config),
      "initial_risks_not_directly_comparable_across_downstream_scales=TRUE",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/evaluate_source_projection_candidates.R")),
      paste0("input_", 2:3, "_manifest_sha256=", vapply(args[2:3],function(root)
        roce_sha256_file(file.path(root,"sha256.txt")),""))), file.path(stage,"metadata.txt"))
    files <- list.files(stage, full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  }, caller = "source candidate heldout evaluation")
  print(reports,row.names=FALSE)
  if (any(reports$failed)) stop("candidate failures retained; no partial selection allowed")
}
if (sys.nframe()==0L) main()

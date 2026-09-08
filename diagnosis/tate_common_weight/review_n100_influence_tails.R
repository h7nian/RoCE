#!/usr/bin/env Rscript

# Saved-data only: reconstruct the primary common-weight TATE pseudo-values
# and site-centered variance, retaining every high-leverage observation.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: review_n100_influence_tails.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- args[[2L]]
  if (file.exists(output)) stop("tail review output already exists")
  source("scripts/slurm/atomic_output.R")
  source("scripts/slurm/result_provenance.R")
  references <- read.csv(file.path(root, "summaries/n100/input_references.csv"))
  stopifnot(identical(references$sim_id, 10001:10100))
  summary_hash <- roce_sha256_file(file.path(root, "summaries/n100/sha256.txt"))
  stopifnot(summary_hash ==
    "bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe")
  cells <- fits <- vector("list", 100L)
  for (i in 1:100) {
    directory <- file.path(root, sprintf("seed_%06d", 10000L+i))
    manifest <- file.path(directory, "sha256.txt")
    stopifnot(roce_sha256_file(manifest) == references$bundle_checksum_manifest_hash[i])
    lines <- readLines(manifest)
    artifact_path <- file.path(directory, "artifacts.rds")
    expected <- substr(lines[substring(lines, 67L) == "artifacts.rds"], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(artifact_path) == expected)
    bundle <- readRDS(artifact_path)
    seed_cells <- seed_fits <- list()
    for (rho in c(0, .5, 1, 1.5, 2, 2.5)) {
      artifact <- bundle$group_result$artifacts[[as.character(rho)]]
      f <- artifact$direct_tate_results$one_round_crossfit
      data <- artifact$data_split
      sizes <- c(t = f$intermediates$sample_sizes$n_t,
                 f$intermediates$sample_sizes$n_source)
      stopifnot(identical(names(sizes), c("t", "s1", "s2")),
                all(sizes == 1000), sum(sizes) == f$N_all,
                all(is.finite(f$all_phi_tau)))
      phi <- counts <- fold_ids <- lapply(sizes, function(n) numeric(n))
      for (k in seq_along(f$intermediates$fold_info)) {
        info <- f$intermediates$fold_info[[k]]
        eta <- f$fold_weights[k, ]
        target <- (1-sum(eta)) * (info$varphi_ot + info$fold_target_estimate)
        for (j in 1:2) {
          target <- target + eta[j] * (info$zeta_components[[j]] + info$mu_pred_ts[j])
          index <- info$source_idx[[j]]
          site <- names(sizes)[j+1L]
          phi[[site]][index] <- f$N_all/sizes[[site]] * eta[j] *
            (info$xi_components[[j]] + info$delta_ts[j])
          counts[[site]] <- counts[[site]] + tabulate(index, nbins = sizes[[site]])
          fold_ids[[site]][index] <- k
        }
        index <- info$target_idx
        phi$t[index] <- f$N_all/sizes[["t"]] * target
        counts$t <- counts$t + tabulate(index, nbins = sizes[["t"]])
        fold_ids$t[index] <- k
      }
      stopifnot(all(unlist(counts) == 1), all(unlist(fold_ids) %in% 1:5))
      assembled <- unlist(phi, use.names = FALSE)
      phi_error <- max(abs(assembled-f$all_phi_tau))
      estimate_error <- abs(mean(assembled)-f$estimate)
      ss <- vapply(phi, function(x) sum((x-mean(x))^2), numeric(1L))
      variance <- sum(ss)/f$N_all^2
      variance_error <- abs(variance-f$se^2)
      stopifnot(phi_error < 1e-10, estimate_error < 1e-10,
                variance_error < 1e-12, variance > 0)
      seed_fits[[length(seed_fits)+1L]] <- data.frame(
        sim_id = 10000L+i, rho, estimate = f$estimate, se = f$se,
        phi_error, estimate_error, variance_error)
      for (site in names(phi)) {
        centered <- phi[[site]]-mean(phi[[site]])
        order <- order(abs(centered), decreasing = TRUE)
        index <- order[1L]
        seed_cells[[length(seed_cells)+1L]] <- data.frame(
          sim_id = 10000L+i, rho, site, site_n = sizes[[site]],
          site_mean_phi = mean(phi[[site]]), site_ss = ss[[site]],
          site_variance = ss[[site]]/f$N_all^2, total_variance = variance,
          row = index, fold = fold_ids[[site]][index],
          A = data[[site]]$A[index], Y = data[[site]]$Y[index],
          raw_phi = phi[[site]][index], centered_phi = centered[index],
          top1_site_ss_fraction = if (ss[[site]] == 0) NA_real_ else
            centered[index]^2/ss[[site]],
          top5_site_ss_fraction = if (ss[[site]] == 0) NA_real_ else
            sum(centered[order[1:5]]^2)/ss[[site]],
          top1_total_variance_fraction = centered[index]^2/sum(ss))
      }
    }
    cells[[i]] <- do.call(rbind, seed_cells)
    fits[[i]] <- do.call(rbind, seed_fits)
    rm(bundle, artifact, f, data)
    if (i %% 10L == 0L) cat("tail-reconstructed seeds:", i, "\n")
  }
  cells <- do.call(rbind, cells)
  fits <- do.call(rbind, fits)
  stopifnot(nrow(cells) == 1800L, nrow(fits) == 600L)
  extremes <- cells[order(abs(cells$centered_phi), decreasing = TRUE), ]
  roce_write_atomic_directory(output, function(stage) {
    write.csv(cells, file.path(stage, "site_tail_cells.csv"), row.names = FALSE)
    write.csv(fits, file.path(stage, "reconstruction_checks.csv"), row.names = FALSE)
    write.csv(head(extremes, 30L), file.path(stage, "largest_tail_cells.csv"), row.names = FALSE)
    writeLines(c("n100_primary_tate_tail_review=numerical_reconstruction_passed",
      "seed_bundles=100", "fits=600", "site_cells=1800",
      paste0("summary_checksum_manifest_fingerprint=", summary_hash),
      paste0("review_script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/review_n100_influence_tails.R")),
      paste0("max_phi_error=", max(fits$phi_error)),
      paste0("max_variance_error=", max(fits$variance_error)),
      paste0("zero_site_ss_cells=", sum(cells$site_ss == 0)),
      "high_leverage_observations_removed=0", "inference_validated=FALSE"),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "saved primary influence tail review")
  print(head(extremes, 10L), row.names = FALSE)
}

if (sys.nframe() == 0L) main()

#!/usr/bin/env Rscript
# Read-only conversion of archived calibrated RoCE fits. Oracle conditional
# moments are diagnostic inputs, not feasible confidence certificates.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L,
          startsWith(arguments[2L], "/scratch.global/zhan9381/FACE-HD/"),
          !dir.exists(arguments[2L]))
root <- arguments[1L]
output <- arguments[2L]
dir.create(output, recursive = TRUE)
manifest <- read.csv(file.path(root, "manifest.csv"))
transform <- rbind(c(0, 0, 1, -1), c(1, 0, -1, 0), c(0, 1, 0, -1))
records <- vector("list", nrow(manifest))
for (index in seq_len(nrow(manifest))) {
  directory <- file.path(root, "cases", manifest$task_id[index])
  receipt <- readRDS(file.path(directory, "COMPLETE.rds"))
  actual_hashes <- tools::md5sum(file.path(directory, names(receipt$hashes)))
  stopifnot(identical(unname(actual_hashes), unname(receipt$hashes)))
  saved <- readRDS(file.path(directory, "population_packets.rds"))
  variant <- saved$variants$common_outcome
  packet <- variant$packet
  population <- variant$population
  source_observed <- do.call(rbind, lapply(packet$source_messages, function(message) message$raw_sum / message$n))
  source_expected <- do.call(rbind, lapply(population$source_moments, `[[`, "mean"))
  observed_common <- sweep(packet$means[-1L, , drop = FALSE], 2L, source_observed[1L, ], "-")[1L, ]
  expected_common <- population$means[2L, ] - source_expected[1L, ]
  stopifnot(max(abs(sweep(packet$means[-1L, , drop = FALSE] - source_observed, 2L, observed_common, "-"))) < 1e-10,
            max(abs(packet$training_shifts)) < 1e-10)
  target_observed <- drop(transform %*% c(packet$means[1L, ], observed_common))
  target_expected <- drop(transform %*% c(population$means[1L, ], expected_common))
  target_truth <- c(expected_common[1L] - expected_common[2L], population$truth - expected_common)
  source_bias <- sweep(source_expected, 2L, target_truth[2:3], "-")
  target_size <- packet$evaluation_sizes[["t"]]
  source_sizes <- unname(packet$evaluation_sizes[-1L])
  source_covariance <- lapply(population$source_moments, `[[`, "covariance")
  empirical_covariance <- Map(function(value, size) value * size / (size - 1),
                              packet$source_private_covariances, source_sizes)
  valid <- variant$valid[-1L, , drop = FALSE]
  source_budget <- vapply(1:2, function(arm) max(abs(source_bias[valid[, arm], arm])), numeric(1L))
  records[[index]] <- list(case_id=manifest$task_id[index], setting=as.list(manifest[index, ]),
    source=unname(source_observed), source_expected=unname(source_expected),
    variance=unname(do.call(rbind, lapply(source_covariance, diag))),
    empirical_variance=unname(do.call(rbind, lapply(empirical_covariance, diag))),
    target=unname(target_observed), target_expected=unname(target_expected),
    target_truth=unname(target_truth), truth=unname(population$truth[1L] - population$truth[2L]),
    target_covariance=unname(transform %*% population$target_covariance %*% t(transform)),
    target_empirical_covariance=unname(transform %*% packet$common_covariance %*% t(transform) * target_size / (target_size - 1)),
    source_bias_budget=unname(source_budget), target_bias_budget=unname(abs(target_expected-target_truth)),
    valid_minimum=unname(colSums(valid)), evaluation_sizes=c(target_size, source_sizes),
    source_packet_md5=unname(tools::md5sum(file.path(directory, "population_packets.rds"))))
}
jsonlite::write_json(records, file.path(output, "packets.json"), auto_unbox=TRUE, digits=17, pretty=TRUE)
writeLines(c("90 archived fits converted without modification.",
             "Conditional population moments and valid-source identities are oracle diagnostic inputs."),
           file.path(output, "SCOPE.txt"))
cat("EXPORTED_FITTED_PACKETS", length(records), "\n")

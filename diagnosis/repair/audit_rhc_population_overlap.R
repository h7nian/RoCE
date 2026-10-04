#!/usr/bin/env Rscript
# Necessary population moment-feasibility checks under bounded merged weights.
# A violated coordinate is a certificate; passing these checks is not a proof
# that all moments can be matched jointly or by the working weight model.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
template_path <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
template <- readRDS(template_path)
features <- template$raw_features
witnesses <- cbind(intercept = 1, features)
binary <- vapply(seq_len(ncol(features)), function(j) all(features[, j] %in% c(0, 1)), logical(1L))
complements <- 1 - features[, binary, drop = FALSE]
colnames(complements) <- paste0("not_", colnames(complements))
witnesses <- cbind(witnesses, complements)
groups <- list(income_reference = grep("^income", colnames(features)),
  race_reference = grep("^race", colnames(features)),
  diagnosis_reference = grep("^cat1", colnames(features)),
  cancer_reference = match(c("caNo", "caYes"), colnames(features)))
for (name in names(groups)) {
  columns <- groups[[name]]
  stopifnot(length(columns) > 0L, !anyNA(columns))
  reference <- 1 - rowSums(features[, columns, drop = FALSE])
  stopifnot(all(reference %in% c(0, 1)))
  witnesses <- cbind(witnesses, reference)
  colnames(witnesses)[ncol(witnesses)] <- name
}
positive <- pmax(witnesses, 0); negative <- pmin(witnesses, 0)
records <- list()
for (scenario in names(template$populations)) {
  population <- template$populations[[scenario]]
  target <- population$sites$t
  for (site_name in setdiff(names(population$sites), "t")) for (arm in 0:1) {
    site <- population$sites[[site_name]]
    joint <- site$probability * if (arm == 1L) site$propensity else 1 - site$propensity
    for (stage in c("initial", "outcome_derivative")) {
      h <- if (stage == "initial") rep(1, length(joint)) else
        target$outcome[, arm + 1L] * (1 - target$outcome[, arm + 1L])
      target_moment <- colSums(witnesses * (target$probability * h))
      positive_moment <- colSums(positive * (joint * h))
      negative_moment <- colSums(negative * (joint * h))
      for (radius in c(2, 3, 5)) {
        maximum <- exp(radius) * positive_moment + exp(-radius) * negative_moment
        minimum <- exp(-radius) * positive_moment + exp(radius) * negative_moment
        gap <- pmax(target_moment - maximum, minimum - target_moment, 0)
        tolerance <- 1e-10 * pmax(1, abs(target_moment), abs(maximum), abs(minimum))
        records[[length(records) + 1L]] <- data.frame(scenario = scenario, site = site_name,
          arm = arm, stage = stage, radius = radius, witness = colnames(witnesses),
          target_moment = target_moment, feasible_minimum = minimum, feasible_maximum = maximum,
          gap = gap, violated = gap > tolerance)
      }
    }
  }
}
result <- do.call(rbind, records)
# True joint weights are strictly inside radii3/5 in W scenarios. Every
# coordinate must pass for any positive derivative function in those cells.
valid_weight_cells <- result$scenario != "O_case_mix" & result$radius >= 3
stopifnot(!any(result$violated[valid_weight_cells]))
dir.create(output, recursive = TRUE)
write.csv(result, file.path(output, "coordinate_moments.csv"), row.names = FALSE)
write.csv(result[result$violated, ], file.path(output, "infeasibility_certificates.csv"), row.names = FALSE)
summary <- aggregate(violated ~ scenario + site + arm + stage + radius, data = result, FUN = sum)
names(summary)[names(summary) == "violated"] <- "violated_coordinates"
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
writeLines(c("Population coordinate feasibility audit for the frozen hypothetical laws, not observed RHC identification.",
  "A violated coordinate proves its calibration moment cannot be matched by any weight in [exp(-M),exp(M)].",
  "No violations do not establish joint feasibility or a full nuisance theorem.",
  "In outcome_derivative checks, h is the target true outcome derivative: correct for the OR-correct O_case_mix limit; an arbitrary positive witness function in the W scenarios."),
  file.path(output, "scope.txt"))
print(summary[summary$violated_coordinates > 0L, ])

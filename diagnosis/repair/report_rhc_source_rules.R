#!/usr/bin/env Rscript
# Render checked, paired results without re-fitting or selecting a favorable split.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4L)
study <- normalizePath(arguments[1L], mustWork = TRUE)
comparison <- normalizePath(arguments[2L], mustWork = TRUE)
weights <- normalizePath(arguments[3L], mustWork = TRUE)
output <- arguments[4L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"),
          !dir.exists(output), file.exists(file.path(comparison, "CHECKS_PASSED")),
          file.exists(file.path(weights, "CHECKS_PASSED")))
manifest <- read.csv(file.path(study, "manifest.csv"), stringsAsFactors = FALSE)
configuration <- jsonlite::fromJSON(file.path(study, "configuration.json"))
results <- read.csv(file.path(comparison, "comparisons.csv"), stringsAsFactors = FALSE)
site_variance <- read.csv(file.path(weights, "site_variance.csv"), stringsAsFactors = FALSE)
weight_summary <- read.csv(file.path(weights, "source_weight_summary.csv"), stringsAsFactors = FALSE)
balance <- read.csv(file.path(weights, "weighted_balance.csv"), stringsAsFactors = FALSE)
results$profile <- paste0("RHC_", results$variant, "_fold", results$fold_seed)
stopifnot(setequal(results$profile, manifest$config),
          setequal(unique(site_variance$profile), manifest$config),
          setequal(unique(balance$profile), manifest$config),
          !anyDuplicated(results$profile))
for (index in seq_len(nrow(results))) {
  name <- results$profile[index]
  sites <- site_variance[site_variance$profile == name, ]
  source_weights <- weight_summary[weight_summary$profile == name, ]
  selected_balance <- balance[balance$profile == name, ]
  gaps <- selected_balance$weighted_minus_target_smd
  undefined <- is.na(gaps)
  if (any(undefined)) stopifnot(all(selected_balance$target_sd[undefined] == 0),
                                length(selected_balance$weighted_minus_target_mean) == length(gaps))
  stopifnot(nrow(sites) == 4L, nrow(source_weights) == 6L,
            abs(sum(sites$variance) - results$se[index]^2) < 1e-12,
            any(!undefined), all(is.finite(gaps[!undefined])))
  results$largest_record_variance_share[index] <- max(sites$top1_total_variance_share)
  results$min_source_ess_fraction[index] <- min(source_weights$effective_n/source_weights$n)
  results$max_source_weight[index] <- max(source_weights$max_weight)
  results$balance_max_abs_smd[index] <- max(abs(gaps[!undefined]))
  results$balance_mean_abs_smd[index] <- mean(abs(gaps[!undefined]))
  results$undefined_balance_smd[index] <- sum(undefined)
  results$max_constant_feature_mean_gap[index] <- if (any(undefined))
    max(abs(selected_balance$weighted_minus_target_mean[undefined])) else 0
}
variants <- unique(results$variant)
paired <- lapply(variants, function(variant) {
  selected <- results[results$variant == variant, ]
  reference <- results[results$variant == "min", ]
  reference <- reference[match(selected$fold_seed, reference$fold_seed), ]
  stopifnot(identical(selected$fold_seed, reference$fold_seed))
  data.frame(variant = variant, partitions = nrow(selected),
    estimate_min_pp = 100 * min(selected$estimate), estimate_max_pp = 100 * max(selected$estimate),
    se_min_pp = 100 * min(selected$se), se_max_pp = 100 * max(selected$se),
    median_se_ratio = median(selected$se/reference$se),
    lower_se_partitions = sum(selected$se < reference$se),
    worst_record_variance_share = max(selected$largest_record_variance_share),
    lower_record_share_partitions = sum(selected$largest_record_variance_share < reference$largest_record_variance_share),
    min_source_ess_fraction_min = min(selected$min_source_ess_fraction),
    min_source_ess_fraction_max = max(selected$min_source_ess_fraction),
    worst_abs_balance_smd = max(selected$balance_max_abs_smd),
    median_abs_balance_smd = median(selected$balance_mean_abs_smd),
    undefined_smd_count = max(selected$undefined_balance_smd),
    max_constant_feature_mean_gap = max(selected$max_constant_feature_mean_gap))
})
summary <- do.call(rbind, paired)
dir.create(output, recursive = TRUE)
write.csv(results, file.path(output, "paired_diagnostics.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
labels <- c(min = "Original min", final_weight_1se = "Final weight 1se",
  final_outcome_1se = "Final OR 1se", initial_weight_1se = "Initial weight 1se",
  source_all_1se = "All source stages 1se", global_1se = "Global 1se")
draw <- function() {
  par(mfrow = c(2, 2), mar = c(4, 10, 3, 1), oma = c(2, 0, 1, 0))
  seeds <- sort(unique(results$fold_seed))
  colors <- grDevices::hcl.colors(length(seeds), "Dark 3")
  for (field in c("estimate", "se", "largest_record_variance_share", "min_source_ess_fraction")) {
    values <- results[[field]] * 100
    titles <- c(estimate = "TATE estimate", se = "Reported standard error",
      largest_record_variance_share = "Largest single-record variance share", min_source_ess_fraction = "Smallest source-arm ESS fraction")
    units <- if (field %in% c("estimate", "se")) "Percentage points" else if (field == "largest_record_variance_share") "Percent of total reported variance" else "Percent of the corresponding source-arm sample"
    plot(range(values), c(.5, length(variants) + .5), type = "n", yaxt = "n",
      xlab = units, ylab = "", main = titles[[field]])
    axis(2, at = seq_along(variants), labels = labels[variants], las = 1, cex.axis = .8)
    abline(h = seq_along(variants), col = "grey90")
    for (variant_index in seq_along(variants)) {
      selected <- which(results$variant == variants[variant_index])
      positions <- variant_index + seq(-.16, .16, length.out = length(selected))
      segments(min(values[selected]), variant_index, max(values[selected]), variant_index, col = "grey65")
      points(values[selected], positions, pch = 19, cex = .65,
        col = colors[match(results$fold_seed[selected], seeds)])
    }
  }
  mtext(sprintf("%d matched fold partitions of one cohort; points are partitions, lines show their range. No coverage claim.",
    length(seeds)), outer = TRUE, side = 1, cex = .8)
}
pdf(file.path(output, "source_regularization.pdf"), width = 13, height = 8)
draw()
dev.off()
# Render the vector figure without an X11 display on compute nodes.
converter <- Sys.which("gs")
stopifnot(nzchar(converter))
status <- system2(converter, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r120",
  shQuote(paste0("-sOutputFile=", file.path(output, "source_regularization.png"))),
  shQuote(file.path(output, "source_regularization.pdf"))), env = "LD_LIBRARY_PATH=")
stopifnot(status == 0L, file.exists(file.path(output, "source_regularization.png")))
table_rows <- vapply(seq_len(nrow(summary)), function(i) {
  row <- summary[i, ]
  sprintf("| %s | %.2f–%.2f | %.2f–%.2f | %.1f%% | %.1f–%.1f |",
    row$variant, row$estimate_min_pp, row$estimate_max_pp, row$se_min_pp, row$se_max_pp,
    100 * row$worst_record_variance_share, 100 * row$min_source_ess_fraction_min, 100 * row$min_source_ess_fraction_max)
}, character(1L))
writeLines(c("# RHC source 正则化配对比较", "",
  sprintf("完整纳入 %d 项拟合、%d 个分折；同一批 5,039 人，target 为 %s，特征配置为 %s，另三个 insurance groups 为 sources。", nrow(results), length(unique(results$fold_seed)), configuration$target_site, configuration$covariate_profile),
  "Two-layer、one-round、joint TATE、cutoff=2、100 个 lambda、半径 5 和求解精度保持一致。", "",
  "下表为跨分折的实际范围，估计和 SE 的单位为百分点；病例方差占比取所有 site、所有分折中的最大值。", "",
  "| 规则 | TATE 范围 | SE 范围 | 最大单病例方差占比 | 各分折最小 source-arm ESS 比例范围（%） |",
  "|---|---:|---:|---:|---:|", table_rows, "",
  "`source_all_1se` 仅改变 source 初始权重、最终权重及最终 OR 的 CV 选取规则；target 及初始 OR 消息保持 min。",
  "`global_1se` 还改变 target 及初始 OR 消息。`1se` 是 nuisance CV 的 one-standard-error 规则，不是修改 aggregation cutoff。", "",
  "[可导出的四面板图](source_regularization.pdf)；[逐分折诊断](paired_diagnostics.csv)；[汇总](summary.csv)。", "",
  "检查包括相同 cohort/fold、实际选取规则、固定阶段数值相等以及报告方差与 influence 贡献重构一致。",
  "普通协变量 SMD 也完整保留：更强正则化可能牺牲平衡；导数加权 calibration 并不要求这里的普通 SMD 精确为零。", "",
  "target 内无变异的特征没有可定义的 target-SD 标准化差值；这些 SMD 标为 NA，并另存原始加权均值差，不用极小分母生成虚大数值。汇总 SMD 仅针对可定义的部分。", "",
  "这些是一个观测队列的不同分折，不能估计 coverage/RMSE，也不能证明偏差更小、消除了混杂或迁移假设成立。",
  "保留所有分折及原始结果；不按显著性筛选，不在本报告中更换默认设置或稿件结果。"),
  file.path(output, "report_zh.md"))
print(summary)

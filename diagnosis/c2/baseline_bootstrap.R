#!/usr/bin/env Rscript
# Multiplier (wild) bootstrap variance for the POOLING baselines (sample_size, IVW).
#
# Diagnosis-only. The analytic variance for sample_size / inverse_variance adds a
# DerSimonian-Laird random-effects tau^2 (estimators_helpers.R calculate_dl_heterogeneity),
# which conflates covariate-shift BIAS with sampling variance and inflates the SE
# (empirically sd(bias)/se ~ 0.6 => CIs too wide, masking bias). This probe recomputes
# their SE by a multiplier bootstrap on the per-site AIPW influence functions
# (no nuisance re-fitting), which yields the honest fixed-effects sampling SE.
# Per seed it records the estimate, truth, analytic SE (with tau^2) and bootstrap SE,
# and coverage under each. RoCE and target_only are NOT bootstrapped (handled elsewhere).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("expected 7 args: tag config n_total K p seed n_folds", call. = FALSE)
tag <- args[[1L]]; config <- args[[2L]]
n_total <- as.integer(args[[3L]]); K_sites <- as.integer(args[[4L]])
p <- as.integer(args[[5L]]); seed <- as.integer(args[[6L]]); n_folds <- as.integer(args[[7L]])
B <- as.integer(Sys.getenv("C2_BOOT_B", "1000"))

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) library(RoCE)
  else if (requireNamespace("devtools", quietly = FALSE)) devtools::load_all(".")
  else stop("Neither installed RoCE nor devtools available.", call. = FALSE)
})
ff <- function(name) { ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) stop(sprintf("missing %s", name), call. = FALSE)
  get(name, envir = ns, inherits = FALSE) }
for (fn in c("generate_simulation_data","split_data_by_site",
             "estimate_sample_size_weighted","estimate_inverse_variance_weighted")) {
  if (!exists(fn, mode="function")) assign(fn, ff(fn))
}
fit_all_sites <- ff(".fit_site_aipw_all_sites")
`%||%` <- function(x, y) if (is.null(x) || !is.finite(suppressWarnings(as.numeric(x)[1]))) y else x

out_root <- Sys.getenv("C2_BOOT_OUTPUT_ROOT", file.path("diagnosis","c2","baseline_bootstrap"))
out_dir <- file.path(out_root, "results"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf("%s_%s_n%d_K%d_p%d_seed%d", tag, config, n_total, K_sites, p, seed))

shift_strength <- as.numeric(Sys.getenv("C2_BOOT_SHIFT", "0.5"))
set.seed(seed)
data <- generate_simulation_data(n_total=n_total, K=K_sites, p=p, config=config,
  estimand_type="superpopulation", site_allocation="model", transform_type="mild",
  outcome_type="binary", heterogeneity_type="none", shift_strength=shift_strength,
  dgp_type="roce", warn_ignored=FALSE)
split <- split_data_by_site(data)
truth <- as.numeric(data$mu1_true)

# --- analytic baselines (carry the DL tau^2) ---
ss  <- estimate_sample_size_weighted(split, family="binomial", use_crossfit=TRUE, n_folds=n_folds, A_val=1L)
ivw <- estimate_inverse_variance_weighted(split, family="binomial", use_crossfit=TRUE, n_folds=n_folds, A_val=1L)

# --- per-site AIPW influence functions for the multiplier bootstrap ---
site_fits <- fit_all_sites(data_split=split, family="binomial", use_rcal=FALSE,
                           use_crossfit=TRUE, n_folds=n_folds, A_val=1L)
sites <- names(site_fits)
mu_j   <- sapply(sites, function(s) as.numeric(site_fits[[s]]$estimate))
n_j    <- sapply(sites, function(s) as.numeric(site_fits[[s]]$n))
ifs    <- lapply(sites, function(s) as.numeric(site_fits[[s]]$varphi_ot))  # centered IF, length n_j
names(ifs) <- sites
V_j_n  <- sapply(sites, function(s) as.numeric(site_fits[[s]]$variance))   # = V_ot/n_j (site sampling var)

# fixed weights
w_ss  <- n_j / sum(n_j)
prec  <- 1 / pmax(V_j_n, 1e-12); w_ivw <- prec / sum(prec)

# multiplier bootstrap: mu_j^(b) = mu_j + mean_i( omega_{j,i} * ifs_j_i ), omega ~ Rademacher (mean0,var1)
boot_ss <- numeric(B); boot_ivw <- numeric(B)
for (b in seq_len(B)) {
  mu_b <- vapply(sites, function(s) {
    v <- ifs[[s]]; nb <- length(v)
    om <- sample(c(-1, 1), nb, replace = TRUE)
    mu_j[[s]] + mean(om * v)
  }, numeric(1))
  boot_ss[b]  <- sum(w_ss  * mu_b)
  boot_ivw[b] <- sum(w_ivw * mu_b)
}

# se_fe = EXACT analytic fixed-effects IF variance (Sigma w_j^2 V_j/n_j), already computed
# in components$var_fixed_effects. This is the proposed exact fix; se_boot is the cross-check.
mk <- function(method, est, se_analytic, se_fe, se_boot) data.frame(
  tag=tag, config=config, n_total=n_total, K=K_sites, p=p, shift=shift_strength, seed=seed, method=method,
  estimate=as.numeric(est), truth=truth, bias=as.numeric(est)-truth,
  se_analytic=as.numeric(se_analytic), se_fe=as.numeric(se_fe), se_boot=as.numeric(se_boot),
  covered_analytic = abs(as.numeric(est)-truth) <= 1.96*as.numeric(se_analytic),
  covered_fe       = abs(as.numeric(est)-truth) <= 1.96*as.numeric(se_fe),
  covered_boot     = abs(as.numeric(est)-truth) <= 1.96*as.numeric(se_boot),
  stringsAsFactors=FALSE)

df <- rbind(
  mk("sample_size",      ss$estimate,  ss$se  %||% sqrt(ss$variance),
     sqrt(ss$components$var_fixed_effects),  sd(boot_ss)),
  mk("inverse_variance", ivw$estimate, ivw$se %||% sqrt(ivw$variance),
     sqrt(ivw$components$var_fixed_effects), sd(boot_ivw))
)
write.csv(df, paste0(prefix, "_boot.csv"), row.names = FALSE)
cat(sprintf("[boot] %s config=%s n=%d B=%d\n", tag, config, n_total, B))
print(df[, c("method","bias","se_analytic","se_fe","se_boot","covered_analytic","covered_fe","covered_boot")], row.names=FALSE)
cat("[boot] DONE\n")

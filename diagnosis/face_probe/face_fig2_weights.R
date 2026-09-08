#!/usr/bin/env Rscript
# FACE.pdf Figure-2 analog: the RoCE adaptive ensemble weights eta-hat for each
# source site across the deviation gradient. As one source (s1) is made
# increasingly non-informative (rho up), its eta-hat should shrink toward 0
# (negative-transfer control) while the informative source(s) keep weight.
# Extracts run_crossfit()$average_weights via load_all + fork; M_tau_inference=5.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(parallel))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "200"))
KK      <- as.integer(Sys.getenv("ROCE_PROBE_K", "2"))
pp      <- as.integer(Sys.getenv("ROCE_PROBE_P", "10"))
n_site  <- 1000L
rho_grid <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5)

message(sprintf("[fig2] eta-hat weights | p=%d C1 K=%d n=%d binary | nsims=%d cores=%d | M_tau=5",
                pp, KK, n_site*(KK+1L), n_sims, n_cores))

one <- function(args) {
  rho <- args$rho; sim <- args$sim
  set.seed(sim)
  d <- generate_face_data(n_total = n_site*(KK+1L), K = KK, p = pp, config = "C1",
                          outcome_type = "binary", estimand_type = "superpopulation",
                          ate_deviation = rho, n_deviated_sites = if (rho > 0) 1L else 0L)
  ds <- split_data_by_site(d)
  r <- tryCatch(run_crossfit(ds, n_folds = 5L, communication_mode = "one_round",
                             lambda_selection = "cv", verbose = FALSE,
                             family = "binomial", M_tau_inference = 5),
                error = function(e) NULL)
  if (is.null(r) || is.null(r$weights)) return(NULL)
  w <- as.numeric(r$weights)                    # average eta-hat, ordered s1, s2, ... (s1 = deviated)
  data.frame(rho = rho, sim = sim,
             source = paste0("s", seq_along(w)), eta = w, stringsAsFactors = FALSE)
}
grid <- do.call(rbind, lapply(rho_grid, function(r) data.frame(rho = r, sim = seq_len(n_sims))))
res <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))), one, mc.cores = n_cores))

agg <- aggregate(eta ~ rho + source, data = res, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
out <- data.frame(rho = agg$rho, source = agg$source,
                  eta_mean = agg$eta[, "mean"], eta_se = agg$eta[, "se"])
out$role <- ifelse(out$source == "s1", "deviated (s1)", "informative")
out$K <- KK
out_csv <- sprintf("diagnosis/face_probe/validation/fig2_weights_K%d.csv", KK)
write.csv(out, out_csv, row.names = FALSE)
message(sprintf("[fig2] wrote %s (%d rows)", out_csv, nrow(out)))

message("\n  rho  | " , paste(sprintf("eta(%s)", unique(out$source)), collapse = "  "))
for (r in rho_grid) {
  s <- out[out$rho == r, ]
  message(sprintf("  %.1f  | %s", r, paste(sprintf("%.3f", s$eta_mean[order(s$source)]), collapse = "    ")))
}
message("[fig2] DONE")

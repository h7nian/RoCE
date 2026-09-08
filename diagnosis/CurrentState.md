# Current state of RoCE (Robust Federated Causal Estimation)

> **Single source of truth for what is true RIGHT NOW.**
> This file is rewritten (not appended) on every successful iteration.
> History lives in `HISTORY.md`.
>
> Last updated: 2026-09-07
> Last updating HISTORY entry: [#0002](HISTORY.md#0002)

---

## 1. Method default (1-line summary per knob)

- Estimand: target average treatment effect (TATE) `tau_t = mu^1_t - mu^0_t`.
- Primary estimator: `run_tate_crossfit(communication_mode = "one_round")`, one
  common source-weight vector per outer fold (FACE eq. 9 structure).
- Nuisances: exponential-tilting site/treatment model and GLM outcome model,
  both fitted on the working basis `phi(X)`; calibrated losses (Tan 2020);
  two-level cross-fitting, 5 outer folds; `nlambda_init = 100`, `lambda` rule
  `min`; truncation `M_tau = M_tau_inference = 5`.
- Aggregation weights: `N_all * Var(eta) + sum_j (lambda t_j - 1)_+ |eta_j|`,
  `lambda = AGG_WALD_LAMBDA = 1` (cutoff c = 1), soft-threshold solve.
- Variance: site-centered pseudo-value variance with weights treated as fixed
  (`se`). **Known deficit: omits the weight-learning term** (see §3).
- Simulation DGP: FACE-style skew-normal covariates, binary outcome,
  `rho = ate_deviation` log-odds treatment shift on source `s1`'s treated arm.
  Working bases still differ across configurations (C2/C3 drop the quadratic
  block from one basis) — **being replaced** (see §3, #0002).

> Authoritative spec: `docs/main.tex` + `docs/supplemental.tex` (Rule 13,
> amendment A1). This section is a 1-line-per-knob restatement.

## 2. Headline numbers (multi-seed verified)

C1, p = 100, K = 2, 1000 observations per site, 100 independent seeds
(v19 package; `results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`):

| rho | RoCE RMSE | target-only RMSE | RoCE coverage | analytic SE / empirical SD | RoCE bias |
|---|---|---|---|---|---|
| 0 | 0.0250 | 0.0348 | 0.93 | 1.00 | -0.012 |
| 0.5 | 0.0259 | 0.0348 | 0.87 | 0.93 | +0.011 |
| 1 | 0.0321 | 0.0348 | 0.82 | 0.77 | +0.011 |
| 1.5 | 0.0332 | 0.0348 | 0.83 | 0.75 | -0.000 |
| 2 | 0.0317 | 0.0348 | 0.90 | 0.83 | -0.007 |
| 2.5 | 0.0301 | 0.0348 | 0.93 | 0.92 | -0.010 |

No C2/C3, K=4/8, or RHC results exist for the current TATE estimator.

## 3. Active issues (open iterations)

1. **#0001 Weight-layer influence for the soft-threshold rule** — Stage 1.
   Symptom: SE omits Var(eta_hat) term. Status: PROPOSED. Blocker: none.
2. **#0002 Common working basis + X-dagger misspecification pilot** — Stage 1.
   Symptom: C2/C3 bases break the common-phi orthogonality assumption of
   `docs/main.tex:222`. Status: PROPOSED. Blocker: none.
3. #0003 Truncation alignment of the tilting calibrated loss (`density_ratio_weight`
   ignores `M_tau`; IF uses `exp(-T_M)`) — Stage 1. Not started.
4. #0004 Inner-fold weight learning on calibrated nuisances (paper: "same
   secondary-fold scheme"; code uses initial plug-ins) — Stage 1, cost gate.
5. #0005 [MIGRATION] of #0001-#0004 into `R/` + `src/`, tests, R CMD check,
   substitute Rule 7a review, Rule 24 audit.
6. #0006 Method-row schema freeze (target_only, roce A, roce B, roce armwise,
   SS, IVW, federated_dr, pooled_dr), shared-shift scenario S, manifest
   rebuild (500 x C1-C3 x K=2,4,8 + S), staged gates n = 10 / 50 / 100.
7. #0007 Production runs, aggregation, figures/tables, paper updates
   (substitute Rule 7c review before numbers leave).
8. #0008 RHC with the frozen package; #0009 cleanup (rename `direct_tate`
   -> `tate`, retire root `main.R`/`realdata.R` legacy pipeline, README,
   archive `diagnosis/tate_common_weight`).

## 4. Parked for later (not v1 blockers)

- Worker warning capture (`warning_capture_complete = FALSE`): structured
  capture candidate exists under `results/.../condition_capture_candidate.*`;
  not deployed.
- DR CV path tail skipping (about 17% of low-lambda candidates invalid per
  replication): add a "selected lambda at valid-path boundary" diagnostic.
- `AGG_WALD_LAMBDA` read from `ROCE_AGG_WALD_LAMBDA` at load time
  (reproducibility hazard); make it an explicit argument.
- One-round vs two-round communication comparison for the common-weight TATE.

## 5. Standing rules (or link to skill RULES.md)

Project-specific amendments (recorded in HISTORY #0001, Rule 27(f)):

## 6. How to resume work (resume protocol)

## 7. Project-specific resources

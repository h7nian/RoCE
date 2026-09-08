# Current state of RoCE (Robust Federated Causal Estimation)

> **Single source of truth for what is true RIGHT NOW.**
> This file is rewritten (not appended) on every successful iteration.
> History lives in `HISTORY.md`.
>
> Last updated: 2026-09-07
> Last updating HISTORY entry: [#0002](HISTORY.md#0002) (DECIDED-PASS)

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
  Working bases still differ across configurations in the package (C2/C3 drop
  the quadratic block); the replacement (common basis, X-dagger
  misspecification at strength 0.75) is validated (#0002) and awaits
  migration (#0006).

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

Target-only coverage 0.92. Mechanism (n=100 decomposition): at rho=1 source
`s1` keeps mean fold weight 0.15 with Wald statistic 1-3, adding +0.023 bias;
fold-weight SD 0.05-0.14 times discrepancy 0.13-0.20 is the missing SD.

**#0001 → #0005 (CLOSED 2026-09-08): delta-method weight-layer SE for the
production soft-threshold rule, now the package's reported `se`, replayed on
the same 100 seeds** — SE/SD 1.03 / 1.03 / 0.99 / 0.96 / 0.97 / 1.02, coverage
0.93 / 0.91 / 0.88 / 0.94 / 0.93 / 0.96 for rho = 0 … 2.5; points unchanged;
package equals the prototype to 6e-17.

**#0007 (CLOSED 2026-09-08): smooth quadratic-bias rule B on the same 100
seeds** — bias −0.011 / +0.006 / +0.007 / +0.004 / +0.002 / +0.001, coverage
0.93 / 0.92 / 0.93 / 0.94 / 0.94 / 0.93; reported alongside rule A as
`<method>_ate_quadratic_bias`.

Common-basis DGP (#0006) single-seed C1–C4 refits (seed 20001, K = 2):
weight-layer SE 0.0220 / 0.0223 / 0.0219 / 0.0214 vs fixed-weight SE 0.0216 /
0.0218 / 0.0216 / 0.0214. No multi-seed C2/C3, K=4/8, or RHC results exist for
the current TATE estimator yet (production #0010).

## 3. Active issues (open iterations)

1. **#0001 Weight-layer influence for the soft-threshold rule** — Stage 1
   DECIDED-PASS (HISTORY #0001). Next: Stage 2 migration (#0005a).
2. **#0002 Common working basis + X-dagger misspecification pilot** — Stage 1
   DECIDED-PASS: omega* = 0.75; C1 unchanged (2e-4). Next: Stage 2 (#0006).
3. **#0003 Truncation alignment of the tilting calibrated loss** — Stage 1
   IN-FLIGHT (candidate library; score check passed, C1 identity refit and
   full tests pending).
4. #0004 Inner-fold weight learning on calibrated nuisances: proposal is to
   align the Supplement text with the implemented out-of-two-fold initial
   plug-ins and run a bounded calibrated-inner sensitivity on C1/K2 instead
   of a x4 nuisance cost in production — awaiting the user's answer.
5. #0005 [MIGRATION] weight-layer variance into `R/` — CLOSED (DECIDED-PASS
   2026-09-08; full suite clean in run 18573008, replay 18573009 exact).
6. #0006 [MIGRATION] common-basis DGP with `misspecification_strength = 0.75`
   — CLOSED (DECIDED-PASS 2026-09-08; four-configuration reproduction
   18573010 exact).
7. #0007 [MIGRATION] smooth quadratic-bias weight rule (`screening_rule =
   "quadratic_bias"`, simulation row `<method>_ate_quadratic_bias` gated by
   `include_quadratic_bias_rule = TRUE`, quadratic weight layer) — CLOSED
   (DECIDED-PASS 2026-09-08; replay 18573011 exact).
8. **#0008 [MIGRATION] truncation-aligned tilting loss into `src/`** — IN-FLIGHT
   (patch applied to main 2026-09-08; #0003 (e)/(f) job 18576245 running; commits
   with the production-library gate of #0010).
9. #0009 [MIGRATION] schema freeze + shared-shift scenario — CLOSED
   (DECIDED-PASS 2026-09-08, merged into main): `.tate_production_method_rows()`,
   `deviation_mechanism = "both_arms"` (both-arm reuse refit), manifests with a
   `deviation_mechanism` column and a `shared_shift` family (C1, K = 2/4/8),
   staged gates pre-registered in HISTORY #0009 §4.
10. **#0010 Production runs** (HISTORY #0010: production root
    `results/direct_tate_mc500_b5000/production_20260908_v1/`, gates A1–A5,
    per-setting grouped submissions with `ROCE_SETTING`, rungs
    1/5/10/25/50/100/200/300/400/500) — starts after #0008 lands; #0011
    aggregation, figures/tables, paper updates (substitute Rule 7c review
    before numbers leave).
11. #0012 RHC with the frozen package; #0013 cleanup (rename `screening_rule`
    -> `weight_rule`, `direct_tate` -> `tate`, retire root `main.R`/`realdata.R`
    legacy pipeline, README, archive `diagnosis/tate_common_weight`).

## 4. Parked for later (not v1 blockers)

- Continuous outcome under C2/C4 with `misspecification_strength > 0` is
  unvalidated and heavy-tailed (#0006 review); production continuous cells are
  C1-only. Add a continuous C2 sanity check before any such cell is run.

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

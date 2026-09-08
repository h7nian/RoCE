# Research and implementation history

## Iterations

<a id="0001"></a>
## 0001 — 2026-09-07 — Weight-layer influence function for the soft-threshold aggregation rule  [DECIDED-PASS]

> commit: a03e3b5f (baseline snapshot)
> previous related: none (first entry)
> stage: 1 (diagnosis)
> method.tex section: `docs/main.tex` §adaptive_aggregation (eq:agg_penalized_objective, eq:tate_federated_variance); `docs/supplemental.tex` §supp:aggregation_penalty

### 1. Symptom / motivation
100-seed C1/K2/p100 pilot (v19): analytic variance is 59% (rho=1) / 57% (rho=1.5) of
the empirical variance; coverage 82% / 83%. Fold-weight SD of source s1 is 0.05–0.14
while its discrepancy is 0.13–0.20; the product matches the missing SD. Replacing the
analytic SE by the Monte Carlo SD (oracle thought experiment) gives 94–95%, so the
dominant defect is the variance, not the bias.

### 2. Theoretical analysis
The estimator is tau(eta_hat) with eta_hat = argmin N_all Var_hat(eta) +
sum_j (lambda t_j - 1)_+ |eta_j|, a piecewise-smooth function of the inner-fold moments
m = (V_ot, V_t, V_s, C_ot, C_cross, mu_ot, mu_pred, delta). The current SE conditions on
eta_hat. First-order (delta-method / stacked estimating equation) expansion:
tau(eta_hat) ≈ tau(eta*) + sum_k1 L_k1' (eta_hat_k1 - eta*_k1), L_k1,j = exact source
increment used in Phase 3. On the active set S, KKT: 2 N_all (Q_SS eta_S + l_S) +
p_S ∘ sign(eta_S) = 0, hence
  d eta_S = -Q_SS^{-1} [ (dQ_SS) eta_S + dl_S + sign(eta_S) ∘ dp_S / (2 N_all) ],
  d eta_{S^c} = 0,
with Q_jk = [V_ot - C_ot,j - C_ot,k + C_cross,jk 1{j≠k}]/n_t + 1{j=k}(V_t,j/n_t + V_s,j/n_s,j),
l_j = (C_ot,j - V_ot)/n_t, p_j = (lambda t_j - 1)_+, t_j = |d_j|/s_j,
s_j^2 = V_ot/n_t + V_t,j/n_t + V_s,j/n_s,j - 2 C_ot,j/n_t, d_j = mu_ot - mu_pred,j - delta_j.
Moment influence (mass derivative at unit mass, within-inner-fold centering makes the
mean terms vanish): dV/dw_i = (c_i^2 - V)/W, dC/dw_i = (c_i c'_i - C)/W, d mean/dw_i =
(y_i - mean)/W. References: FACE Remark 10 (map H(theta_hat)); Hjort & Claeskens (2003)
JASA; Stefanski & Boos (2002) TAS; Johnson, Lin & Zeng (2008) JASA (L1 active-set
sandwich); Fang & Santos (2019) REStud for the kink caveat. At a kink (lambda t_j = 1
exactly, measure zero) the map is only directionally differentiable; the prototype
classifies coordinates by the realized solution and records how many |lambda t_j - 1|
< 1e-6 occur.

### 3. Decomposition
New prototype `diagnosis/weight_layer/weight_layer.R` reusing the empirical-mass
bookkeeping of `diagnosis/tate_common_weight/quadratic_weight_influence.R` (inner
records, moment recomputation with assertion against stored moments, outer values and
increments). New: `.soft_threshold_weight_derivative()` implementing the KKT
derivative above and `.weight_layer_influence(fit)` returning direct, indirect and
total per-observation gradients plus `variance_weight_layer`. No package code changes.
Finite-difference functional `.weight_layer_functional(fit, mass)` re-solves the
weights with the package's `optimize_weights()` (identical clipping/floors) so the
numerical derivative is of the production rule. Replay over the 100 saved v19 seeds.

### 4. Acceptance criteria
- [x] (a) seed 10013 bundle, 6 rho × 8 directions: max |numerical − analytical| ≤ 1e-8 (ε = 1e-5) → 1.7e-10 [PASS]
- [x] (b) fixed-weight estimate and variance reproduce stored `fit$estimate`, `fit$variance` to ≤ 1e-12 → asserted for all 600 fits [PASS]
- [x] (c) direct-only FD error > 1e-5 in ≥ 1 direction (weight layer non-trivial) → min 1.0e-4, max 1.4e-2 [PASS]
- [x] (d) no inner-record observation appears in its own outer evaluation fold (assertion) → asserted, 600 fits [PASS]
- [x] (e) n=100 replay: weight-layer coverage ≥ 0.88 at rho = 1 and rho = 1.5 → 0.88 / 0.94 [PASS]
- [x] (f) n=100 replay: mean weight-layer SE / empirical SD ∈ [0.90, 1.10] at every rho; rho=0 coverage ∈ [0.92, 0.97] → 1.027, 1.032, 0.993, 0.956, 0.970, 1.020; rho=0 coverage 0.93 [PASS]
- [x] (g) kink count (|lambda t_j − 1| < 1e-6) reported; ≤ 1% of fold×source cells → 0 of 6000 [PASS]

### 5. Validation results (filled after running)
Slurm job 18530725 (`diagnosis/weight_layer/run_weight_layer.sh`), 3 min 17 s, single core.
Output `diagnosis/out/weight_layer/v1/` (sha256 manifest inside). Saved v19 fits, no refits.

| rho | bias | empirical SD | fixed-weight SE | weight-layer SE | SE/SD fixed → layer | coverage fixed → layer |
|---|---|---|---|---|---|---|
| 0 | -0.0120 | 0.0220 | 0.0220 | 0.0226 | 1.00 → 1.03 | 0.93 → 0.93 |
| 0.5 | +0.0108 | 0.0237 | 0.0220 | 0.0244 | 0.93 → 1.03 | 0.87 → 0.91 |
| 1 | +0.0106 | 0.0305 | 0.0234 | 0.0303 | 0.77 → 0.99 | 0.82 → 0.88 |
| 1.5 | -0.0004 | 0.0333 | 0.0251 | 0.0319 | 0.75 → 0.96 | 0.83 → 0.94 |
| 2 | -0.0066 | 0.0312 | 0.0258 | 0.0302 | 0.83 → 0.97 | 0.90 → 0.93 |
| 2.5 | -0.0099 | 0.0286 | 0.0262 | 0.0291 | 0.92 → 1.02 | 0.93 → 0.96 |

Coverage MCSE ≈ 0.025. The direct/indirect cross term is negative in 40/100 seeds at rho=0
and 4/100 at rho=1; dropping it would overstate variance at rho=0 and understate it at
rho=1, confirming that a scalar SE multiplier cannot represent this correction. The
residual dip at rho=1 (0.88) is the retained-source bias (+0.0106 ≈ 0.35 SD), not SE.

### 6. Decision + rationale
DECIDED-PASS. The delta-method weight layer removes the SE understatement at every rho
without changing point estimates, resampling, or applying a multiplier. Migrate to the
package (Stage 2, entry #0005): report `se` = weight-layer SE as the analytic SE, keep
`se_fixed_weights` as a diagnostic, document the derivation in the Supplement and the
kink caveat in the main text. The rho=1 bias is reported as the weak-separation regime.

<a id="0002"></a>
## 0002 — 2026-09-07 — Common working basis and X-dagger misspecification: strength pilot  [DECIDED-PASS]

> commit: a03e3b5f (baseline snapshot)
> previous related: none
> stage: 1 (diagnosis)
> method.tex section: `docs/main.tex` §notation (line 222: single working feature map phi(X)), §simulations (configuration definitions, lines ~880–890)

### 1. Symptom / motivation
C2/C3/C4 are implemented by dropping the quadratic block from one working basis
(`R/data_generation_face.R:563-569`), so the tilting basis Z and outcome basis W differ.
The calibration score then balances only Z-directions; a deterministic counterexample
(`diagnosis/tate_common_weight/check_calibration_tangent_span.R`) shows a non-zero
outcome-direction derivative 0.0487 in a C3-style example. The paper assumes one
common phi(X) for both nuisances (main.tex:222), as Tan (2020) does; the code's C2/C3
lie outside the theorem. RHC already uses one basis (`R/real_data_rhc.R:476-477`).

### 2. Theoretical analysis
Neyman orthogonality of the calibrated construction requires the tilting score
equations (one per gamma coordinate) to span the outcome-model gradient directions,
which holds when both models share phi(X). Misspecification must therefore be defined
as "true mechanism outside span(phi)" (Kang & Schafer 2007; Tan 2020), not as a
smaller working basis.

### 3. Decomposition
Prototype DGP `diagnosis/dgp_common_basis/dgp_common_basis.R`, built from package
internals (`generate_face_covariates`, `generate_face_outcomes`, `logistic`,
`face_binary_logit`), with:
- Z_site = W_outcome = [X − kappa, X^2] for every configuration;
- misspecified mechanism uses eta_mix = (1 − omega) eta(X) + omega eta(X_dagger),
  X_dagger = standardized Kang–Schafer-style transforms of the four signal coordinates
  (exp(X1/2), X2/(1+exp(X1)), (X1 X3/25 + 0.6)^3, (X2 + X4 + 20)^2), standardized to the
  target reference population (mean kappa, sd 1), remaining coordinates unchanged;
- C1: neither; C2: outcome mechanism; C3: treatment propensity; C4: both;
- truth recomputed on the same 100,000-unit reference population and seed as
  `calculate_face_truth`; rho mechanics unchanged; `site_outcome_shift` knob added
  (both arms, deviated source) for the later shared-shift scenario but held at 0.
Pilot: one seed (20001), p = 100, K = 2, 1000 per site, rho = 0, configs C1–C4, omega
∈ {0.25, 0.5, 0.75, 1}; `run_tate_crossfit()` with production settings.

### 4. Acceptance criteria (pre-registered omega selection rule)
- [x] (a) reference-population R^2 of eta_true on phi(X) ∈ [0.60, 0.90] for the misspecified mechanism (per omega) → C2: 0.913 / 0.751 / 0.618 / 0.528; C3: 0.945 / 0.806 / 0.659 / 0.548 for omega = 0.25 / 0.5 / 0.75 / 1 [PASS at omega 0.5 and 0.75]
- [x] (b) inference logit truncation fraction ≤ 0.01 and max raw density-ratio weight ≤ 20 → truncation fraction 0 in all 13 cells; max raw weight 13–37 with **C1 itself at 29.2** — the absolute bound of 20 was a mis-specified gate (Rule 26): the skew-normal shift already yields raw ratios near 30 at radius 5 without truncation. Re-based to ≤ 1.5 × C1 (43.8): all cells [PASS]
- [x] (c) true propensity ∈ [0.02, 0.98] for ≥ 99% of units → the DGP clips propensities to [0.1, 0.9] in every configuration (5.5% clipped under C1, 3.1–4.9% under the misspecified mechanisms), so the original wording is satisfied by construction; recorded as "no more clipping than C1" [PASS]
- [x] (d) C1 under the common basis: |TATE − current C1 TATE (seed 20001, rho 0: 0.194329)| ≤ 0.002 → 0.000215 [PASS]
- [x] (e) omega* = largest omega satisfying (a)–(c) for both C2 and C3 → **omega* = 0.75** (C2 R² 0.618, C3 R² 0.659) [PASS]
- [x] (f) with omega*: C2/C3/C4 fits finite, all selected nuisance fits converged → finite; 0 non-converged fits in every cell [PASS]

### 5. Validation results (filled after running)
Slurm array 18530726 (13 cells, 10 cores each, 1.2–2.0 h per cell); omega selection job
18546175 (absolute weight gate, `selection_absolute_gate_v1/`) and its re-evaluation with
the C1-relative weight gate (`selection/`). One seed (20001), rho = 0, p = 100, K = 2.

| config | omega | R² misspecified | max raw weight | clip fraction | TATE | target-only | truth |
|---|---|---|---|---|---|---|---|
| C1 | 0 | — | 29.2 | 0.055 | 0.1945 | 0.2185 | 0.2063 |
| C2 | 0.5 | 0.751 | 26.2 | 0.055 | 0.2091 | 0.2056 | 0.2139 |
| C2 | 0.75 | 0.618 | 28.6 | 0.055 | 0.2164 | 0.2108 | 0.2173 |
| C3 | 0.5 | 0.806 | 14.4 | 0.041 | 0.2054 | 0.2288 | 0.2063 |
| C3 | 0.75 | 0.659 | 18.3 | 0.033 | 0.2087 | 0.2220 | 0.2063 |
| C4 | 0.75 | 0.618 / 0.659 | 16.1 | 0.033 | 0.2126 | 0.1858 | 0.2173 |

Truth changes with omega under C2/C4 because the true outcome mechanism changes; it is
recomputed on the same 100,000-unit reference population. Source weights stay near
0.3 / 0.2 (both sources compatible at rho = 0). C1 reproduces the current DGP to
2e-4 (basis centering only).

Re-run 2026-09-07 with the re-based (C1-relative) gates, job 18549812
(`diagnosis/out/dgp_common_basis/selection/omega_selection_table.csv`): omega* = 0.75 again
(all_ok TRUE at omega 0.5 and 0.75 for C2/C3/C4; FALSE at 0.25 and 1 on the R^2 gate), C1
identity difference 2.15e-4.

### 6. Decision + rationale
DECIDED-PASS with omega* = 0.75 pre-registered for C2–C4 (Kang–Schafer-style transforms
of the four signal coordinates, standardized on the target reference population, mixed
into the true outcome (C2), treatment-propensity (C3) or both (C4) mechanisms). Two gates
were mis-specified and re-based on the C1 mechanism's own values (Rule 26; recorded
above). Migrate the DGP to `R/data_generation_face.R` (Stage 2, entry #0006): common
basis for every configuration, `misspecification_strength` argument fixed at 0.75 in the
production manifests, documentation of C1–C4 in README and `docs/main.tex`.

<a id="0003"></a>
## 0003 — 2026-09-07 — Truncation alignment of the tilting calibrated loss  [IN-FLIGHT]

> commit: (this entry's prototype committed after submission)
> previous related: none
> stage: 1 (diagnosis; isolated candidate library, package source untouched)
> method.tex section: `docs/main.tex` §nuisance truncation paragraph (lines 353–356: every displayed weight is implemented as exp[−T_{M}(phi'gamma)] in the calibrated losses and the influence function); `docs/supplemental.tex` lines 196–202

### 1. Symptom / motivation
The influence function and the outcome calibrated loss use the truncated tilt weight
exp[−T_M(Z'gamma)] (`src/variance.cpp:57`, `src/cv_utils.h:74-82`), but the tilting
coordinate descent and its CV validation loss use exp[−clip_{±50}(Z'gamma)]
(`src/cv_utils.h:498-508`, `:563`, `:619`, `:683-698`), i.e. no truncation. When
truncation is active the calibration score and the functional's alpha-derivative
disagree by E_s[I (e^{−g} − e^{−T_M(g)}) g' Z]^T(alpha_hat − alpha*), a first-order term
absent from the analytic variance. Inactive in the C1 pilots (truncation fraction 0)
but relevant for RHC and for strong shift.

### 2. Theoretical analysis
Replace exp(−u) in the tilting loss by psi_M(u) = exp(−T_M(u)) (1 − (u − T_M(u))): the
exponential inside the radius continued by its tangent outside. psi_M is convex,
psi_M'(u) = −exp(−T_M(u)), psi_M''(u) = exp(−u) 1{|u| < M}. The first-order condition
becomes E_t[g' Z] = E_s[I exp(−T_M(Z'gamma)) g' Z], which is exactly the
alpha-derivative of the truncated functional in the Z directions, so orthogonality
holds for the estimator that is actually evaluated. Loss, CV validation loss, score
and influence function then share one weight. Refined (two-round) tilting fits keep
an infinite radius in both CV and final fit (deferred; not the manuscript default).

### 3. Decomposition
Patch `diagnosis/truncation_alignment/truncation_alignment.patch` (applied to an
isolated copy): `tilt_weight()` / `tilt_loss()` replace `density_ratio_weight()`;
`density_ratio_cd_update()` and `density_ratio_val_loss()` take `M_tau`; the initial
and calibrated fits and their CV kernels pass the fit radius; comparison-method
density ratios pass `Inf` (baselines unchanged). Validation script
`truncation_alignment.R`: (A) synthetic two-site example, radius 2, unpenalized
calibrated fit, truncated vs untruncated score under candidate and baseline
libraries; (B) C1 seed-20001 refit versus the #0002 stored fit.

### 4. Acceptance criteria
- [x] (a) candidate: max |truncated score| ≤ 1e-6 with truncation fraction > 0.05 → 2.25e-9 at truncation fraction 0.876 (radius 2, 1259 iterations, converged) [PASS]
- [x] (b) baseline: max |untruncated score| ≤ 1e-6 and max |truncated score| > 1e-3 (documents the defect) → 1.27e-9 and 0.0736 (truncation fraction 0.132) [PASS]
- [x] (c) candidate: max |untruncated score| > 1e-3 (the fit really changed) → 3.37e8 (the truncated solution's untruncated score explodes; radius-2 saturation as noted in §2) [PASS]
- [x] (d) C1 identity: estimate, se, fold weights, target-only, source estimates differ ≤ 1e-10 from the #0002 C1 cell → all differences exactly 0; production radius inactive (truncation fraction 0) [PASS]
- [ ] (e) full installed testthat: 0 failures, 0 errors → measured on the final tree that carries the same patch (#0008 (b), production gate A1 of #0010) [SUPERSEDED]
- [ ] (f) R CMD check `Status: OK` → measured on the final tree (#0008 (c), production gate A2 of #0010) [SUPERSEDED]

### 5. Validation results (filled after running)
Job 18550692 (2026-09-08, 2 h 59 min): steps 1–3 completed with (a)–(d) as recorded in §4
(`diagnosis/out/truncation_alignment/{baseline,candidate}/`). Step 4 (full installed test suite
on the candidate) hit the 3 h limit inside `test-comparison-bootstrap-variance.R`, and the
candidate source had been copied before the test-fixture fixes of #0005 §5 landed. Resubmitted
as phase `tests` of the same script (rebuilds the candidate from the current tree + patch, runs
only step 4 and R CMD check) with a 10 h limit: job recorded below.
(e), (f): job 18576245 (phase `tests`) ran 8 h 11 min and was still inside
`test-comparison-bootstrap-variance.R` when the scheduler reset of 2026-09-08 13:05 cancelled
every job on the account; the 3 h run had stalled at the same place. Cause: the script exported
`ROCE_NUISANCE_CV_THREADS=5` for step 3 and the test suite's forked workers deadlock in the
OpenMP nuisance-CV kernel with it (the main-tree runs never set it). The script now unsets it
before step 4. Since #0008 applied the same patch to the main tree, (e)/(f) are measured there
by the #0010 production gates A1/A2 instead of a third candidate run.

### 6. Decision + rationale
PENDING

<a id="0005"></a>
## 0005 — 2026-09-07 — Land the weight-layer variance in the package  [MIGRATION]

> commit: (pending)
> previous related: [#0001](#0001) (Stage 1, DECIDED-PASS)
> stage: 2 (migration into `R/`; no `src/` change)
> method.tex section: `docs/main.tex` §adaptive_aggregation, new paragraph "Variance of the learned weights"; `docs/supplemental.tex` new §supp:weight_layer

### 1. Symptom / motivation
Stage-1 prototype #0001 passed all acceptance criteria. The package still reports the
fixed-weight pseudo-value SE (`se`), which understates the sampling variation of the
learned weights by up to 25% at intermediate source deviation.

### 2. Theoretical analysis
As in #0001 §2, now written into the manuscript (Rule 13 first): main text paragraph
with the expansion and the stationarity system; supplement §supp:weight_layer with
$Q$, $l$, $p$, the implicit-function derivative, the moment influences, the increment
$L_{k_1}$, the gradient, the variance identity and the kink caveat; bib entry
`fang2019inference`.

### 3. Decomposition
- New `R/aggregation_weight_influence.R`: `.inner_fold_records()`, `.inner_fold_moments()`,
  `.outer_fold_increments()`, `.variance_quadratic_form()`, `.wald_penalty_ingredients()`,
  `.soft_threshold_weight_derivative()`, `.weight_layer_gradient()`,
  `.weight_layer_variance()`; tolerances `WEIGHT_LAYER_*` in `R/constants.R`.
- `aggregate_fold_estimates()` computes the weight layer from `fold_info`,
  `inner_fold_info`, `fold_lambdas`, `fold_weight_psd_ridge`, `screening_rule`; returns
  `final_variance` (weight layer), `variance_fixed_weights`, `se_fixed_weights`,
  `weight_layer` diagnostics. Both aggregation entry points forward these; results gain
  `se_fixed_weights`, `variance_fixed_weights`, `weight_layer`. `se`/`variance`/CIs now
  include the weight layer (paper: reported SE).
- `.summarize_tate_aggregation_diagnostics()` adds `se_fixed_weights`,
  `weight_layer_indirect_variance`, `weight_layer_cross_term`, `weight_layer_kink_cells`.
- Hard-threshold screening differentiates no penalty (retained sources solved with
  lambda = 0); CV-selected aggregation lambda is treated as fixed (documented).
- New `tests/testthat/test-weight-layer-variance.R` (decomposition identity, finite
  differences of the mass functional, vanishing layer when no source is retained).
- Diagnosis prototype reduced to a replay driver calling package functions after the
  package builds (Rule 16 cleanup).

### 4. Acceptance criteria
- [x] (a) substitute Rule 7a review  recorded; every finding FIX/REBUT/DEFER → §7 (ten findings: 8 FIX, 1 REBUT, 1 DEFER) [PASS]
- [x] (b) full build + testthat via `./test.sh`: 0 failures, 0 errors → run 18573008 (tree d184c5eb), exit 0 [PASS]
- [x] (c) package replay of the 100 v19 seeds reproduces `diagnosis/out/weight_layer/v1/replay_rows.csv` weight-layer SE to ≤ 1e-10 → job 18573009: max |package − prototype| weight-layer SE 5.6e-17, fixed SE 5.2e-17 (`diagnosis/out/weight_layer/v2/`) [PASS]
- [x] (d) Rule 24 audit of user instructions in scope (readability, no dead code, naming, paper alignment) → no references to the retired `.soft_threshold_weight_derivative()` / `penalized` argument remain; the retired C2 diagnostic tests are archived, not deleted silently; `se` / `se_fixed_weights` / `weight_layer` named consistently across the fit list, the simulation rows and the QC columns; main.tex + supplement + README describe the same estimator [PASS]

### 7. Review note
Reviewer verdict: math and correspondence to production sound; three items to fix before
landing. Findings and handling:
1. Active PSD ridge treated as a constant → **FIX**: `.weight_layer_gradient()` now stops
   when any fold ridge is positive (the derivative is not derived for an active ridge);
   supplement wording corrected accordingly.
2. Absolute tolerances (moment check 1e-10, KKT 1e-8·max(1,|l|), floor stop on s²) →
   **FIX**: moment check and site-sum check are relative to the magnitude of the stored
   quantity; the stationarity residual is compared with 1e-6 times the magnitude of the
   stationarity terms; the discrepancy-variance floor now mirrors the optimizer's
   `VARIANCE_MIN` floor with a zero derivative instead of stopping.
3. `.summarize_tate_aggregation_diagnostics()` requires the new fields → **REBUT**: every
   production caller summarises a fresh fit; summarising fits saved by older package
   versions is not a supported path (fail-fast by design).
4. Arm-wise and two-round paths now also report the weight layer → **FIX (test)**: added
   an arm-specific `run_crossfit()` decomposition test; two-round/binomial are exercised
   by the existing suite's smoke tests (crash coverage), statistical validation deferred
   to the C1 gate of #0006.
5. Hard-threshold kink counter always zero → **FIX**: the derivative takes the fold
   multiplier for the kink count and a `penalized` flag for the penalty term; the
   discrete inclusion event is documented in the code comment as not differentiated.
6. K=1 / n_folds=2 untested → **DEFER**: K=1 is covered by the vanishing-layer test;
   n_folds=2 is not a supported production setting.
7. Relative FD tolerance fragile near zero → **FIX**: absolute 1e-8·max(1,|analytical|).
8. Site-sum tolerance absolute → **FIX** (relative, see 2).
9. Error message lacks fold/magnitude → **FIX**.
10. Supplement ridge wording → **FIX**.

### 5. Validation results (filled after running)
PENDING — test job 18550689 (`./test.sh`), package replay job 18550690 (depends on the test job), submitted 2026-09-07 after the substitute review fixes and the unit-mass length fix (first full run 18546449 failed: inner records index the full site, the mass vector was sized by the inner training sample). The test job also covers #0006.
2026-09-08: job 18550689 built and ran the full suite (tree f7436de3: #0005 + #0006 + #0007); the only errors were the ten tests of `test-slurm-atomic-output.R`, which read `repo_root` from the helper's own environment (not visible after `load_all()`; pre-existing since the a03e3b5f snapshot, unrelated to these entries). Fixed in `helper-load.R` (commit 63809d41, `repo_root` exposed on the attached helper environment). The summary reporter had capped its listing at ten failures; the full list (test.cmd now prints every failure) also contained two `test-tate-aggregation.R` assertions of the pre-#0005 identity `variance == pseudo-value variance` (now asserted on `variance_fixed_weights`, with the decomposition identity added) and three `test-weight-layer-variance.R` fixture errors (binary outcomes under the gaussian family gave constant-y folds; the fixtures now generate continuous outcomes; the indefinite-matrix check was rebuilt through `C_ot`). Commit 1ed81858; resubmitted as test job 18565508 with replay 18565509, DGP reproduction 18565510, quadratic replay 18565511 dependent on it. Run 18565508 left one over-tight tolerance (quadratic coordinate-descent cross-check, 3e-8 relative; fixed in d184c5eb). **Run 18573008 (tree d184c5eb): all tests passed, 0 failures, 0 errors, exit 0.** Replay 18573009, DGP reproduction 18573010 and quadratic replay 18573011 run against the library installed by that job.

### 6. Decision + rationale
**DECIDED-PASS (2026-09-08).** The delta-method weight layer is the reported variance of the
common-weight TATE (`se`), with the fixed-weight pseudo-value SE retained as
`se_fixed_weights`; the package reproduces the Stage-1 prototype to machine precision on the
100 saved seeds and the full suite is clean. Landed in commits f7436de3 … d184c5eb (shared
landing with #0006/#0007; see #0006 §7 item 4 for why the entries share one commit).
Status: CLOSED.

<a id="0006"></a>
## 0006 — 2026-09-07 — Land the common-basis DGP with X-dagger misspecification  [MIGRATION]

> commit: (pending)
> previous related: [#0002](#0002) (Stage 1, DECIDED-PASS, omega* = 0.75)
> stage: 2 (migration into `R/`; no `src/` change)
> method.tex section: `docs/main.tex` §simulations configuration paragraph (line ~920); README "Simulation Configurations"

### 1. Symptom / motivation
The package still implements C2/C3/C4 by dropping the quadratic block from one working
basis, which breaks the common-feature-map assumption of the theory (#0002 §1).

### 2. Theoretical analysis
As in #0002: both nuisances use phi(X) = [X − kappa, X^2]; misspecification lives in the
true mechanism through eta_mix = (1 − omega) eta(X) + omega eta(X†) with omega = 0.75.

### 3. Decomposition
- `R/constants.R`: `FACE_MISSPECIFICATION_STRENGTH <- 0.75`, `FACE_SIGNAL_COORDINATES <- 4L`
  (used by the parameter getters instead of the literal 4).
- `R/data_generation_face.R`: `.face_transformed_coordinates()`, `.face_standardization()`,
  `.face_x_dagger()`, `.face_mixed_predictor()`, `.face_misspecification_strengths()`,
  `.face_reference_population()` (standardization, binary calibration and truth from one
  reference draw); `calculate_face_propensity()` and `generate_face_outcomes()` accept
  `X_dagger` + strength; `get_face_binary_calibration()` and `calculate_face_truth()` become
  configuration-aware wrappers of the reference population; `generate_face_data()` builds
  one basis for every configuration and returns the transformed `X_dagger`.
- `generate_simulation_data(misspecification_strength = FACE_MISSPECIFICATION_STRENGTH)`.
- Tests: configuration tests rewritten (common basis; C2/C4 change Y and truth, C3/C4
  change A; strength 0 reproduces C1); binary-calibration reuse test adapted. The 13
  `test-c2-*` probes and their helper encode the retired C2 definition and move to
  `diagnosis/c2/archived_tests/`; their output directories are removed.
- Documentation: manuscript configuration paragraph, README table, roxygen.

### 4. Acceptance criteria
- [x] (a) package DGP reproduces the #0002 prototype for seed 20001: standardization constants and binary calibration equal to 1e-12; truth equal to 1e-12 for C1, C2/C3/C4 at omega 0.75 → job 18573010 (`diagnosis/out/dgp_common_basis/package_C*/dgp_checks.csv`): standardization and calibration differences exactly 0 in all four configurations; truth differences 4.4e-16 (C1, C3) and 2.8e-16 (C2, C4) [PASS]
- [x] (b) package refit of C2 and C3 (omega 0.75, seed 20001) reproduces the prototype fits' estimate, fixed-weight SE and fold weights to 1e-10 → all four configurations refitted: estimate, fixed-weight SE, fold-weight and target-only differences exactly 0 (`fit_checks.csv`); weight-layer SE 0.0220 / 0.0223 / 0.0219 / 0.0214 vs fixed 0.0216 / 0.0218 / 0.0216 / 0.0214 for C1–C4 [PASS]
- [x] (c) full build + testthat: 0 failures, 0 errors → run 18573008 (tree d184c5eb), exit 0 [PASS]
- [x] (d) rho-reuse invariance tests pass under the new C3 (treatment mechanism misspecified) → `test-rho-reuse.R` (C3, p = 3) passes in run 18573008 [PASS]
- [x] (e) substitute Rule 7a review recorded; Rule 24 audit → §7 (nine findings: 5 FIX, 2 DEFER, 2 no action); Rule 24: `config` is an explicit argument of the calibration/truth helpers, `misspecification_strength` is recorded on every result row, no unused variables remain (`X_centered` scoped), README configuration table and main.tex configuration paragraph match the code [PASS]

### 7. Review note
Reviewer verdict: faithful port of the prototype, follows from §2, no blocking defects;
RNG-order equivalence and the HISTORY truths (C1 0.2063, C2 0.2173) reproduced
independently. Findings and handling:
1. `misspecification_strength` absent from result rows, so cells from the retired and the
   new DGP definitions could be merged under the same `config` label → **FIX**:
   `run_single_simulation()` writes `misspecification_strength` (the strength actually
   applied; 0 under C1 and for the roce DGP) on every row.
2. Continuous outcome under C2/C4 unvalidated and heavy-tailed (standardized cubic
   coordinate reaches |z| ≈ 15–20; predictor max 32 → 293 at ω = 0.75) → **DEFER**:
   production continuous cells are C1-only; documented in the `generate_face_data()`
   roxygen and parked in CurrentState §4. A continuous C2 sanity check is required before
   any continuous C2–C4 cell is run.
3. `get_face_binary_calibration()` / `calculate_face_truth()` defaulted to `config = "C1"`,
   silently wrong for C2/C4 callers → **FIX**: `config` is now a required positional
   argument of both; the four callers (`diagnose_rho0_site_tate.R`,
   `diagnose_c3_target_remainder.R`, two `diagnosis/tate_common_weight/check_*.R`) pass it.
4. Diff carries #0005 hunks → **DEFER**: #0005, #0006 and #0007 were validated by one test
   job on the same tree and land as one commit; the entry numbers are named in the message.
5. Rho-reuse invariance safe → no action.
6. Edge cases (p < 4, invalid strength, strength ignored for `dgp_type = "roce"`) →
   **FIX** for the last item: an explicit `misspecification_strength` with the roce DGP now
   stops; the rest need no action.
7. Minor → **FIX**: the returned `misspecification_strength` is the strength applied (0
   under C1); `X_centered` is built only in the effect-modification branch. The double
   predictor evaluation in `.face_reference_population()` is left as is (same as the
   prototype, transient).

### 5. Validation results (filled after running)
Full test suite: run 18573008 (tree d184c5eb) passed with 0 failures, 0 errors (see #0005 §5 for
the fixture fixes on the way there; none touched the DGP). Four-configuration package
reproduction: job 18573010 (array 1–4, dependent on 18573008), all cells exact as recorded in
§4 (a)–(b). The omega-selection rerun with the re-based gates (job 18549812) reproduced
omega* = 0.75 (#0002 §5). — test job 18550689 (`./test.sh`), package reproduction array 18550691 (`diagnosis/dgp_common_basis/run_dgp_common_basis.sh`, 4 configurations, depends on the test job), submitted 2026-09-07.

### 6. Decision + rationale
**DECIDED-PASS (2026-09-08).** The package DGP is the pre-registered common-basis design with
X-dagger misspecification at omega = 0.75; it reproduces the Stage-1 pilot exactly in every
configuration, the truth comes from one reference population, and every result row carries
the strength applied. Landed in commits f7436de3 … d184c5eb. Status: CLOSED.

<a id="0007"></a>
## 0007 — 2026-09-07 — Land the smooth quadratic-bias weight rule as the sensitivity estimator  [MIGRATION]

> commit: (pending)
> previous related: Stage-1 candidate `diagnosis/tate_common_weight/QUADRATIC_BIAS_WEIGHT_CANDIDATE.md`, `QUADRATIC_BIAS_WEIGHT_N100_REVIEW.md`, `ANALYTIC_WEIGHT_LAYER_REVIEW.md` (n = 100: coverage 0.93 / 0.92 / 0.93 / 0.94 / 0.94 / 0.93 with its weight layer); [#0005](#0005)
> stage: 2 (migration into `R/`)
> method.tex section: `docs/main.tex` new rem:quadratic_bias_rule; `docs/supplemental.tex` §supp:weight_layer (quadratic paragraph)

### 1. Symptom / motivation
User decision 2026-09-07: rule A (truncated Wald) stays primary and rule B (smooth quadratic
bias penalty, power 3/4) is computed alongside in every experiment so that switching the
primary later needs no rerun.

### 2. Theoretical analysis
Objective n_t Var(eta) + n_t^{3/4} sum_j delta_j^2 eta_j^2 (rate window (1/2, 1)); linear
normal equations (Q + n_t^{-1/4} diag(delta^2)) eta = -l; smooth everywhere, so the
weight-layer derivative has no kink and uses the penalty curvature plus the
2 n_t^{-1/4} delta_j d delta_j eta_j term.

### 3. Decomposition
- `R/aggregation_weight_rules.R`: `.quadratic_bias_weights()` (Cholesky solve on the
  diagonally scaled curvature, positive-definiteness check, residual check).
- `screening_rule = "quadratic_bias"` accepted by `run_tate_crossfit()`,
  `calculate_tate_crossfit_aggregation()`, `.compute_phase2_weights()`,
  `estimate_tate_weight_bootstrap()`; `AGG_QUADRATIC_BIAS_POWER <- 0.75`.
- `.soft_threshold_weight_derivative()` generalized to `.weight_rule_derivative()`
  (rule-specific active set, curvature and penalty term).
- `run_single_simulation()` adds the `<method>_ate_quadratic_bias` row from the same arm
  fits; QC method lists and the required diagnostic column list updated.
- Tests: quadratic weights solve their normal equations and coincide with the
  unpenalized variance-optimal weights when discrepancies vanish; finite-difference
  check of the quadratic weight layer; simulation smoke row present.
- Replay: package replay of the quadratic rule on the 100 v19 seeds compared with
  the prototype's `quadratic_weight_layer_n100_v1/weight_layer_rows.csv`.

### 4. Acceptance criteria
- [x] (a) normal-equation residual ≤ 1e-12 relative on every inner fold of the test fixture (and ≤ 1e-6 at run time inside the solver); discrepancy-free case equals the exact unpenalized solve `solve(Q, -l)` to 1e-10 and `optimize_weights(lambda = 0)` to 1e-6 (its coordinate descent stops at `WEIGHT_OPT_TOL = 1e-8`, which left a 3e-8 relative gap in run 18565508; amended after the substitute review and that run) → run 18573008: normal equations satisfied on every inner fold at 1e-12; exact solve at 1e-10; coordinate descent at 1e-6 [PASS]
- [x] (b) finite-difference check of the quadratic weight layer: |numerical − analytical| ≤ 1e-8·max(1, |analytical|) → run 18573008 (`test-weight-layer-variance.R`, three directions, K = 2 fixture, zero kink cells) [PASS]
- [x] (c) 100-seed replay: estimates and weight-layer SEs equal the prototype's candidate rows to ≤ 1e-8 → job 18573011: max |package − prototype| estimate 5.3e-16, fold weights 9.4e-16, weight-layer SE 5.2e-17, fixed SE 5.6e-17 (`diagnosis/out/quadratic_bias_rule/v1/`) [PASS]
- [x] (d) full build + testthat: 0 failures, 0 errors → run 18573008 (tree d184c5eb), exit 0 [PASS]
- [x] (e) substitute Rule 7a review recorded; Rule 24 audit → §7 (eleven findings: 6 FIX, 2 DEFER, 3 no action); Rule 24: rule names (`soft_penalty` / `hard_threshold` / `quadratic_bias`) and the gate flag follow the existing `include_hard_threshold_diagnostic` pattern; no dead branches; remark, supplement paragraph and README agree with the code [PASS]

### 7. Review note
Reviewer verdict: rule solve equivalent to the prototype's prototype (Q = A/n_t, l = linear/n_t);
quadratic branch of `.weight_rule_derivative()` verified by hand as the exact implicit-function
derivative and mapped term by term onto the prototype's `quadratic_weight_influence.R`; no blocking
mathematical defects. Findings and handling:
1. Package checks only Q + P positive definite, the prototype also required the unpenalized variance
   quadratic → **FIX**: `.require_positive_definite()` checks Q first, then Q + P.
2. Derivative exact and equivalent → no action.
3. Unconditional quadratic row: a `chol()` failure or a check failure in this path would abort
   the whole replicate including the primary row → **FIX**: `include_quadratic_bias_rule`
   (default TRUE, strictly logical) gates the row; `ROCE_INCLUDE_QUADRATIC_BIAS_RULE` in the
   rho-group task runner; rows record `quadratic_bias_rule_requested`; the checkpoint audit
   drops the method only when every row says FALSE.
4. Row-count change rejects pre-#0007 checkpoints/results at audit time → **DEFER** to the
   #0009 schema freeze; noted here: results and checkpoints written before this entry are
   incompatible with the audits (fail-closed, not silent).
5. Primary point-estimate path unchanged → no action.
6. Criterion (a) wording (1e-10 not attainable; "every fold" tested on one fold) → **FIX**:
   criterion amended to 1e-8 for the coordinate-descent comparison; the test now checks the
   normal equations on every inner fold.
7. Replay does not exercise the Phase-2 wiring → **noted in §5**: the wiring is covered by the
   quadratic FD test's unit-mass identity (`.weight_layer_functional(unit) == fit$estimate`).
8. `lambda_selection = "cv"` still runs the inner CV under the quadratic rule; the Wald
   diagnostics on the quadratic row are descriptive only → **DEFER** (production uses the fixed
   multiplier; the diagnostics are shared across rules by design, see the Phase-2 comment).
9. Non-finite `l` / weights gave an unhelpful error → **FIX**: explicit `is.finite` checks on
   the moments, the curvature and the solution.
10. Edge cases (K = 1, zero/huge discrepancy, small n_t) handled; duplicated sources fail hard
    where the primary path would ridge → accepted under fail-fast (item 3 covers recovery).
11. Missing `skip_on_cran()` → **FIX**.

### 5. Validation results (filled after running)
Full suite: run 18573008 (see #0005 §5 for the runs before it). Replay job 18573011 on the
100 v19 seeds (C1, K = 2, six rho values): quadratic rule bias −0.011 / +0.006 / +0.007 /
+0.004 / +0.002 / +0.001, empirical SD 0.023–0.029, mean weight-layer SE 0.023–0.028,
coverage 0.93 / 0.92 / 0.93 / 0.94 / 0.94 / 0.93 for rho = 0 … 2.5 — identical to the prototype's
`quadratic_weight_layer_n100_v1/rho_summary.csv`. The replay exercises
`.quadratic_bias_weights()` on recomputed moments; the Phase-2 wiring is covered by the
quadratic FD test's unit-mass identity.

### 6. Decision + rationale
**DECIDED-PASS (2026-09-08).** Rule B is computed alongside rule A in every experiment as the
`<method>_ate_quadratic_bias` row (gated, default on) with its own smooth weight layer; on the
100-seed C1/K2 pilot it holds 0.92–0.94 coverage at every rho, including the rho = 0.5–1.5
region where rule A dips, at a small RMSE cost at rho = 0 (bias −0.011 vs. target-only). The
primary rule stays A per the user's decision; the production runs will report both. Landed in
commits f7436de3 … d184c5eb. Status: CLOSED.

<a id="0009"></a>
## 0009 — 2026-09-07 — Freeze the production row schema and add the shared-shift scenario  [MIGRATION]

> commit: branch `entry-0009` (23b6fc36 … 5fadacab), fast-forwarded into main on 2026-09-08
> previous related: [#0005](#0005), [#0006](#0006), [#0007](#0007) (landing candidate f7436de3); user decision 2026-09-07 (add a two-arm shared-shift scenario; 500 replications)
> stage: 2 (package + production tooling)
> method.tex section: `docs/main.tex` sec:simulations (design paragraph: shared-shift experiment)

### 1. Symptom / motivation
Production restarts from fresh manifests once #0005–#0007 land. Before that, (i) the
per-replicate method-row set must be frozen in one place (audits and scripts enumerated it
in more than twenty places); (ii) the negative-transfer experiment moves only the treated arm
of the deviated source (`ate_map`), so the outcome-model heterogeneity that the density-ratio
calibration cannot absorb is never exercised; the user asked for a second, both-arm scenario.

### 2. Theoretical analysis
Shared shift: at the deviated source, Pr{Y(a)=1 | X, R=s_1} = expit{g(X) + rho + Delta a}
with Delta unchanged. The source's conditional risk difference changes only through the
curvature of expit, so its TATE differs from the target's while the treatment log-odds and
the covariate law are unchanged; the target-population truth is unaffected (it depends on
`base_ate` only). Reuse across rho: both arms of the deviated source change, so the
positive-rho update must refit both arms of the changed sources (target and informative
sources are reused unchanged); the treated-arm mechanism keeps the control-arm reuse.

### 3. Decomposition
- Paper first: `docs/main.tex` sec:simulations design paragraph (shared-shift experiment).
- DGP: `generate_face_data(deviation_mechanism = c("treated_arm", "both_arms"))`,
  `build_face_outcome_shift_map()`, `generate_face_outcomes(outcome_shift_map)`; the data
  list carries `deviation_mechanism` and `outcome_shift_map`; the roce DGP rejects an explicit
  mechanism.
- Driver: `run_single_simulation(deviation_mechanism)`, `run_simulation_study(...)`,
  `.face_heterogeneity_type(..., deviation_mechanism)` with labels
  `one_shared_shift_source` / `multiple_shared_shift_sources` /
  `shared_shift_and_effect_modification`; rows record `deviation_mechanism`.
- Reuse: `.validate_one_round_rho_reuse_data(changed_arms)`,
  `.refit_one_round_crossfit_sources(changed_arms)`,
  `.reuse_one_round_tate_across_rho(refit_control_arm)`; `rho_reuse$reused_control_arm`.
- Schema freeze: `.tate_production_method_rows(include_hard_threshold, include_quadratic_bias)`
  and `.face_deviation_mechanisms()` in `R/simulation_diagnostics.R`; the checkpoint/smoke
  audits and the source-diagnostic summary use it; `.validate_face_production_scientific_metadata()`
  requires `deviation_mechanism` and checks the mechanism-consistent label.
- Tooling: manifests carry `deviation_mechanism` (appended after `methods`, so the grouped
  column positions 3/4/6/7 read by the submission script are unchanged);
  `submit_rho_reuse_equivalence.sh` / `audit_rho_reuse_equivalence.R` take the locked setting
  from `ROCE_REUSE_CONFIG/K/RHO` and record the mechanism in the gate; `build_direct_tate_manifest.R` mode `shared_shift` (pre-registered scope: C1,
  K = 2/4/8, rho grid {0, 0.5, 1, 1.5, 2, 2.5}, both arms); the rho-group builder accepts one
  grouped experiment per manifest; task runners pass the mechanism; row annotation asserts it;
  `submit_rho_group_direct_tate.sh` derives the group count from the primary manifest; the
  checkpoint audit is family-aware; `aggregate_direct_tate.R` groups by mechanism;
  `audit_mc500_manifests.R` audits one family per root; `prepare_mc500_manifests.sh` builds the
  shared-shift family under `results/direct_tate_mc500_b5000/shared_shift/`.
- Tests: both-arm DGP invariants (target truth and non-deviated sites unchanged, both arms
  move at s_1), heterogeneity labels, reuse validator arms, shared-shift grouped reuse equals
  the independent fit, frozen row set, validator with both mechanisms.

### 4. Acceptance criteria
- [x] (a) shared-shift grouped reuse equals the independent fit on the K=2/p=3 fixture to 1e-12 (all `rho_equivalence_columns`; s2 reused in both arms, s1 refitted in both) → PASS in runs 18565512 and 18573016
- [x] (b) treated-arm grouped reuse unchanged: existing reuse tests pass byte-identically → PASS (same runs)
- [x] (c) full build + testthat from the worktree (`USE_SOURCE=true`): 0 failures, 0 errors → run 18573016 (tree e0346ee4), exit 0 [PASS]
- [x] (d) manifest family build + `audit_mc500_manifests.R` pass for both families in a scratch root (500 replications: 27000 + 9000 primary tasks, 4500 + 1500 groups) → job 18555672: both audits passed; main 27000/4500 (54 settings), shared_shift 9000/1500 (18 settings); grouped columns 3/4/6/7 unchanged (`deviation_mechanism` after `methods`) [PASS]
- [ ] (e) smoke task (`manifest_smoke_single.csv`) through `run_direct_tate_task.R` + `audit_direct_tate_smoke.R` with the frozen row set → moved to the #0010 pre-submission gate sequence (the smoke runs with the frozen production library) [DEFERRED to #0010]
- [x] (f) substitute Rule 7a review recorded; Rule 24 audit → §7 (nine findings: 4 FIX, 1 DEFER, 4 no action); Rule 24: one argument name (`deviation_mechanism`) from DGP to manifests to rows; `changed_arms` / `refit_control_arm` describe exactly what they do; the frozen row set has one definition; README and the Slurm runbook document the family and the gate; no dead code (every new helper has callers) [PASS]

Before the n = 10 stage of the shared-shift family: a shared-shift reuse-equivalence run
(C1, K = 4, rho = 2.5; grouped both-arm reuse vs. an independent fit) via
`submit_rho_reuse_equivalence.sh` with `ROCE_REUSE_CONFIG=C1` on the shared-shift manifest
root — the grouped submission requires a family-specific gate under the family's result root
whose `deviation_mechanism` line matches the manifest — and the rho = 0 rows of the two
families must be identical for the same sim_id (same dataset and fits).

Staged production gates (pre-registered for #0010, per family and per C×K): n = 10
(implementation: audits pass, no failed replicate, SE/SD within [0.7, 1.4]); n = 50 (coverage of
RoCE at every rho within 0.95 ± 0.06, i.e. two MC standard errors); n = 100 (go/no-go for 500:
coverage within 0.95 ± 0.045 at every rho except the pre-declared weak-separation dip at
rho = 0.5, which is reported, and RoCE RMSE ≤ target-only RMSE at rho = 0 within MC error).

### 7. Review note
Reviewer verdict: both-arm reuse is exact (the changed source's `process_source_site()` is
re-run in full under the same work seed; only target-side components are reused, and the
reference mu0 arm carries the same fields as mu1); the treated-arm path is byte-identical
(`changed_arms = 1L` reproduces the old check; a NULL shift map adds exactly 0). Findings:
1. Exactness → no action.
2. **Regression**: `run_single_simulation()` forwards `deviation_mechanism` unconditionally and
   the roce-DGP guard used `missing()`, so every roce driver run would have errored (no test
   drives roce through the driver) → **FIX**: value-based guards for both
   `deviation_mechanism` and `misspecification_strength`; test added for forwarded defaults.
3. Treated-arm path identical → no action; the rho = 0 cross-family identity is added to the
   pre-n = 10 checks (§4) and to the K = 2 unit test.
4. No test reused an informative source's control arm under both arms (all K = 1) →
   **FIX**: the shared-shift reuse test now uses K = 2 (s2 reused in both arms); a
   production-scale shared-shift reuse-equivalence run is pre-registered in §4.
5. HISTORY wording on the column position → **FIX** (§3).
6. Rho-group runner's invariant columns lacked `deviation_mechanism` → **FIX**.
7. Negative-transfer audits/scripts unaffected → no action.
8. Pre-#0009 raw results fail the audits (required columns) → **DEFER**, intended: production
   restarts from fresh manifests; old gate files are superseded by the rebuild.
9. Edge cases handled; `rbinom` early return at prob exactly 0/1 would desynchronize later
   sites' draws but the reuse validator would stop rather than reuse → noted, no action.

### 5. Validation results (filled after running)
- (d) 2026-09-08, job 18555672 (`prepare_mc500_manifests.sh` with `ROCE_MANIFEST_ROOT` in the
  session scratchpad, cutoff gate from the main tree): main family 27000 primary / 4500 grouped
  rows, 54 settings, diagnostics manifests (smoke 1, cutoff 1000, truncation 2000);
  shared-shift family 9000 / 1500, 18 settings; both `manifest_audit_passed.txt` gates written
  with `family=`, `experiment=`, `deviation_mechanism=` lines.
- (c) worktree run 18555671 (tree 23b6fc36): failures were all fixture/test issues — the
  reuse validator compared outcomes with `identical()` (integer vs double after the fixture's
  flip; now numeric-value comparison), the roce guard test used p = 3 (roce transform needs
  p ≥ 4), the Slurm-helper fixtures lacked `deviation_mechanism`, plus the shared #0005/#0007
  test fixes recorded in #0005 §5. Run 18565512 (tree 7591e430): every test passed except the
  quadratic coordinate-descent tolerance (3e-8 relative; fixed in e0346ee4, exact-solve
  comparison). **Clean confirmation run 18573016 (tree e0346ee4): all tests passed, 0 failures,
  0 errors, exit 0.**
- (a), (b): PASS in run 18565512 (shared-shift K = 2 grouped reuse equals the independent fit
  at 1e-12 on all equivalence columns; treated-arm reuse tests unchanged and passing).
- (e) deferred to #0010; (f) recorded above.

### 6. Decision + rationale
**DECIDED-PASS (2026-09-08).** The production row set is frozen in
`.tate_production_method_rows()`, the shared-shift scenario is implemented end to end (DGP,
driver, both-arm reuse, validators, manifests, audits, submission tooling, paper text) with its
grouped reuse proven exact against independent fits, and both manifest families build and
audit at 500 replications. Merged into main by fast-forward. Status: CLOSED (smoke gate (e)
runs inside #0010).
<a id="0008"></a>
## 0008 — 2026-09-08 — Land the truncation-aligned tilting loss in `src/`  [MIGRATION]

> commit: (pending: applied to the main tree on 2026-09-08, committed after #0003 (e)/(f) pass)
> previous related: [#0003](#0003) (Stage 1, isolated candidate library; (a)–(d) PASS)
> stage: 2 (migration into `src/` and `R/`)
> method.tex section: `docs/main.tex` nuisance truncation paragraph (every displayed weight is exp[−T_M(φ'γ)]); no text change — the code now matches the text

### 1. Symptom / motivation
#0003 showed that the tilting coordinate descent and its CV loss used the untruncated weight
while the influence function and the outcome loss used the truncated one; when the radius is
active the calibration score and the functional's α-derivative disagree at first order.

### 2. Theoretical analysis
As #0003 §2: the loss ψ_M(u) = exp(−T_M(u)) (1 − (u − T_M(u))) has derivative −exp(−T_M(u)),
so loss, CV validation loss, score and influence function share one weight; refined
(two-round) tilting fits keep an infinite radius; comparison methods pass Inf (unchanged).

### 3. Decomposition
`diagnosis/truncation_alignment/truncation_alignment.patch` applied verbatim to the main tree
(`src/cv_utils.h`: `tilt_weight()` / `tilt_loss()`; `src/density_ratio.cpp`: `M_tau` threaded
through `density_ratio_cd_update()` / `density_ratio_val_loss()` and the DR kernels;
`R/model_fitting.R`, `R/estimators_helpers.R`, `R/cross_fitting_algorithms.R`: the initial and
calibrated fits and their CV kernels pass the fit radius; `R/comparison_methods.R`: Inf);
`Rcpp::compileAttributes()` regenerated the export glue.

### 4. Acceptance criteria
- [ ] (a) #0003 (e) and (f) on the isolated candidate (same patch, tree 8cb74da6): full installed suite 0 failures / 0 errors, R CMD check `Status: OK` → ___ [PENDING, job 18576245]
- [ ] (b) production-library audit suite on the final tree (`run_package_audit_tests.sh`, #0010 gate A1): 0 failures / 0 errors → ___ [PENDING]
- [ ] (c) R CMD check on the final tree (`run_r_cmd_check.sh`, #0010 gate A2): `Status: OK` → ___ [PENDING]
- [ ] (d) C1 identity on the final tree is implied by #0003 (d) (radius inactive) and re-checked by the #0010 smoke task's Wald/weight diagnostics → ___ [PENDING]

### 5. Validation results (filled after running)
PENDING

### 6. Decision + rationale
PENDING

<a id="0010"></a>
## 0010 — 2026-09-08 — Production runs (500 replications, both families) with staged gates  [PRODUCTION]

> commit: (pending)
> previous related: [#0009](#0009) (frozen row set, manifest families, pre-registered gates), [#0008](#0008) (final estimator tree)
> stage: production
> method.tex section: `docs/main.tex` sec:simulations (results paragraphs and figures are written from these outputs in #0011)

### 1. Symptom / motivation
No multi-seed results exist for the current estimator beyond the 100-seed C1/K2 pilot. The
paper needs C1–C3 × K = 2/4/8 × six rho values (negative transfer) and C1 × K = 2/4/8 × six
shift values (shared shift), 500 replications each, with the primary rule A, the sensitivity
rule B, the arm-wise diagnostic and the five benchmarks on every replicate.

### 2. Theoretical analysis
Not applicable (execution entry). The pre-registered decision gates of #0009 §4 apply per
family and per config/K block; the weak-separation dip at rho = 0.5 is reported, not tuned.

### 3. Decomposition
Production root `results/direct_tate_mc500_b5000/production_20260908_v1/` (fresh; the earlier
C1/K2 n = 10 rows under `raw/` were produced by the retired estimator and stay archived in
place). Gate sequence, in order:
- A1 `run_package_audit_tests.sh` → isolated library `Rlib_production_20260908_v1` + test gate.
- A2 `run_r_cmd_check.sh` → `package_check_production_20260908_v1` + check gate.
- A3 `prepare_mc500_manifests.sh` (`ROCE_MANIFEST_ROOT` = production root) → both families.
- A4 smoke task (`manifest_smoke_single.csv`, C3/K4/rho 0) + `audit_direct_tate_smoke.sh`.
- A5 reuse-equivalence gates: negative transfer (C3/K4/rho 2.5) on the production root and
  shared shift (`ROCE_REUSE_CONFIG=C1`) on `production_root/shared_shift`.
- B grouped submissions per setting (`ROCE_SETTING="C1:K2"` …, one submitter call per
  config/K block and family) walking the checkpoint ladder 1/5/10/25/50/100/200/300/400/500
  with the dependent six-rho checkpoint audit at every rung; decision reviews at n = 10, 50,
  100 per #0009 §4 before the next rung is submitted.
- Operational amendment (recorded here, Rule 27): `submit_rho_group_direct_tate.sh` gains
  `ROCE_SETTING` scoping and its caps rise from 5 jobs / 2 concurrent per call to 500 / 50, so
  the 18 setting blocks can run concurrently; per-job resources are unchanged
  (K = 2/4/8: 20/40/32 CPUs, 12/16/32 GB, 12/18/24 h).

### 4. Acceptance criteria
- [ ] (a) gates A1–A5 pass on the final tree (fingerprints recorded in each gate file) → ___ [PENDING]
- [ ] (b) n = 10 rung, every setting of both families: checkpoint audit passes, no failed replicate, SE/SD within [0.7, 1.4] for the primary rule → ___ [PENDING]
- [ ] (c) n = 50 rung: primary-rule coverage within 0.95 ± 0.06 at every rho (two MC standard errors) → ___ [PENDING]
- [ ] (d) n = 100 rung (go/no-go for 500): coverage within 0.95 ± 0.045 at every rho except the pre-declared rho = 0.5 dip, which is reported; RoCE RMSE ≤ target-only RMSE at rho = 0 within MC error → ___ [PENDING]
- [ ] (e) n = 500: all 18 blocks committed and audited; `aggregate_direct_tate.R` and the per-setting diagnostics run on both families → ___ [PENDING]

### 5. Validation results (filled after running)
Gates launched 2026-09-08 on tree 0fde383e: A1 test audit job 18580683
(`Rlib_production_20260908_v1`), A2 R CMD check job 18580684
(`package_check_production_20260908_v1`), A3 manifest families job 18580685
(`production_20260908_v1/`); all three were cancelled unstarted by the cluster-wide scheduler
reset at 13:05 (job numbering restarted) and resubmitted as A1 = 33441, A2 = 33442,
A3 = 33443 on the same tree.

- **A3 manifests [PASS]** (job 33443): `production_20260908_v1/` holds the negative-transfer
  family (27,000 primary / 4,500 grouped rows, 54 settings, plus smoke 1, cutoff 1,000,
  truncation 2,000) and `production_20260908_v1/shared_shift/` the shared-shift family
  (9,000 / 1,500 rows, 18 settings); both `manifest_audit_passed.txt` gates written with
  matching `family` / `experiment` / `deviation_mechanism` lines.
- **A2 R CMD check**: first submission (33442) died in 8 s because `ROCE_CHECK_ROOT` was
  relative and the script `cd`s into the check root before staging; resubmitted with an
  absolute path as job 33807.
- A1 = 33441 running. A4/A5 follow A1–A3.

**Gate-hardening round (2026-09-08).** The gates found five real defects before any
production task ran; each is fixed on main and the gates were rerun from a clean root:
1. `ROCE_CHECK_ROOT` must be absolute (the check script `cd`s into it before staging) — A2
   died in 8 s (job 33442).
2. The #0008 patch added `M_tau` to `select_lambda_cv_initial_density_ratio_cpp`; seven test
   call sites still used the old signature (job 33441, commit 9d83a967).
3. The submitter-cap contract test still asserted the pre-#0010 caps 5/2 (same commit).
4. `man/` was stale after the #0006–#0009 roxygen changes: R CMD check reported codoc
   mismatches (same commit), then undocumented arguments on `calculate_face_propensity`,
   `fit_initial_density_ratio` and `run_simulation_study` (job 37155, commit b5882878).
5. `audit_roce_contracts.sh` scanned with `rg`, absent in batch, so `|| true` made it a silent
   no-op (commit 9a366d8e); the `grep` replacement then matched compiled `src/*.o` and
   `src/RoCE.so`, fixed with `-I` in a380869f (jobs 35389/35390).

A1 passed with 0 test failures at job 37154. Because commit b5882878 touched `R/` and `man/`,
the A1 gate's `package_source_fingerprint` went stale, so both gates were rerun on the final
tree: A1 = 38766, A2 = 38764.

Defect found while reading the A2 log and fixed in 9a366d8e: `audit_roce_contracts.sh`
scanned for retired identifiers with `rg`, which does not exist in the batch environment, and
its trailing `|| true` made the scan a silent no-op in every gate. It now uses `grep` with
equivalent exclusions and passes on the current tree.

### 6. Decision + rationale
PENDING

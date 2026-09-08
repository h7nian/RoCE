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
- [ ] (a) candidate: max |truncated score| ≤ 1e-6 with truncation fraction > 0.05 → ___ [PENDING]
- [ ] (b) baseline: max |untruncated score| ≤ 1e-6 and max |truncated score| > 1e-3 (documents the defect) → ___ [PENDING]
- [ ] (c) candidate: max |untruncated score| > 1e-3 (the fit really changed) → ___ [PENDING]
- [ ] (d) C1 identity: estimate, se, fold weights, target-only, source estimates differ ≤ 1e-10 from the #0002 C1 cell → ___ [PENDING]
- [ ] (e) full installed testthat: 0 failures, 0 errors → ___ [PENDING]
- [ ] (f) R CMD check `Status: OK` → ___ [PENDING]

### 5. Validation results (filled after running)
PENDING — Slurm job 18550692 (resubmitted after an internal-function call fix and the #0005 unit-mass fix; the C1 identity refit now uses the #0006 in-package DGP; `diagnosis/truncation_alignment/run_truncation_alignment.sh`), submitted 2026-09-07.

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
- [ ] (a) substitute Rule 7a review  recorded; every finding FIX/REBUT/DEFER → ___ [PENDING]
- [ ] (b) full build + testthat via `./test.sh`: 0 failures, 0 errors → ___ [PENDING]
- [ ] (c) package replay of the 100 v19 seeds reproduces `diagnosis/out/weight_layer/v1/replay_rows.csv` weight-layer SE to ≤ 1e-10 → ___ [PENDING]
- [ ] (d) Rule 24 audit of user instructions in scope (readability, no dead code, naming, paper alignment) → ___ [PENDING]

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

### 6. Decision + rationale
PENDING

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
- [ ] (a) package DGP reproduces the #0002 prototype for seed 20001: standardization constants and binary calibration equal to 1e-12; truth equal to 1e-12 for C1, C2/C3/C4 at omega 0.75 → ___ [PENDING]
- [ ] (b) package refit of C2 and C3 (omega 0.75, seed 20001) reproduces the prototype fits' estimate, fixed-weight SE and fold weights to 1e-10 → ___ [PENDING]
- [ ] (c) full build + testthat: 0 failures, 0 errors → ___ [PENDING]
- [ ] (d) rho-reuse invariance tests pass under the new C3 (treatment mechanism misspecified) → ___ [PENDING]
- [ ] (e) substitute Rule 7a review recorded; Rule 24 audit → ___ [PENDING]

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
PENDING — test job 18550689 (`./test.sh`), package reproduction array 18550691 (`diagnosis/dgp_common_basis/run_dgp_common_basis.sh`, 4 configurations, depends on the test job), submitted 2026-09-07.

### 6. Decision + rationale
PENDING

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
- [ ] (a) normal-equation residual ≤ 1e-12 relative on every inner fold of the test fixture (and ≤ 1e-6 at run time inside the solver); discrepancy-free case equals `optimize_weights(lambda = 0)` to 1e-8 (its coordinate descent stops at `WEIGHT_OPT_TOL = 1e-8`; amended after the substitute review) → ___ [PENDING]
- [ ] (b) finite-difference check of the quadratic weight layer: |numerical − analytical| ≤ 1e-8·max(1, |analytical|) → ___ [PENDING]
- [ ] (c) 100-seed replay: estimates and weight-layer SEs equal the prototype's candidate rows to ≤ 1e-8 → ___ [PENDING]
- [ ] (d) full build + testthat: 0 failures, 0 errors → ___ [PENDING]
- [ ] (e) substitute Rule 7a review recorded; Rule 24 audit → ___ [PENDING]

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
PENDING — the full test job 18550689 was still queued when #0007 landed in the tree, so it
builds and tests #0005 + #0006 + #0007 together; package replay job 18552835
(`diagnosis/quadratic_bias_rule/run_quadratic_bias_rule.sh`, depends on the test job). The
replay exercises `.quadratic_bias_weights()` on recomputed moments; the Phase-2 wiring is
covered by the quadratic FD test's unit-mass identity.

### 6. Decision + rationale
PENDING

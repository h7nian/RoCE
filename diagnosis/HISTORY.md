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
## 0002 — 2026-09-07 — Common working basis and X-dagger misspecification: strength pilot  [IN-FLIGHT]

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
- [ ] (a) reference-population R^2 of eta_true on phi(X) ∈ [0.60, 0.90] for the misspecified mechanism (per omega) → ___ [PENDING]
- [ ] (b) inference logit truncation fraction ≤ 0.01 and max raw density-ratio weight ≤ 20 in every fitted source arm → ___ [PENDING]
- [ ] (c) true propensity ∈ [0.02, 0.98] for ≥ 99% of units → ___ [PENDING]
- [ ] (d) C1 under the common basis: |TATE − current C1 TATE (seed 20001, rho 0: 0.194329)| ≤ 0.002 → ___ [PENDING]
- [ ] (e) omega* = largest omega satisfying (a)–(c) for both C2 and C3; recorded here and in `docs/main.tex` → ___ [PENDING]
- [ ] (f) with omega*: C2/C3/C4 fits finite, all selected nuisance fits converged → ___ [PENDING]

### 5. Validation results (filled after running)
PENDING — Slurm array 18530726 (`diagnosis/dgp_common_basis/run_dgp_common_basis.sh`, 13 cells), submitted 2026-09-07.

### 6. Decision + rationale
PENDING

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
PENDING — Slurm job 18530797 (`diagnosis/truncation_alignment/run_truncation_alignment.sh`, depends on the #0002 C1 cell), submitted 2026-09-07.

### 6. Decision + rationale
PENDING

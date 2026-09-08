# Research and implementation history

## Iterations

<a id="0001"></a>
## 0001 — 2026-09-07 — Weight-layer influence function for the soft-threshold aggregation rule  [IN-FLIGHT]

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
- [ ] (a) seed 10013 bundle, 6 rho × 8 directions: max |numerical − analytical| ≤ 1e-8 (ε = 1e-5) → ___ [PENDING]
- [ ] (b) fixed-weight estimate and variance reproduce stored `fit$estimate`, `fit$variance` to ≤ 1e-12 → ___ [PENDING]
- [ ] (c) direct-only FD error > 1e-5 in ≥ 1 direction (weight layer non-trivial) → ___ [PENDING]
- [ ] (d) no inner-record observation appears in its own outer evaluation fold (assertion) → ___ [PENDING]
- [ ] (e) n=100 replay: weight-layer coverage ≥ 0.88 at rho = 1 and rho = 1.5 → ___ [PENDING]
- [ ] (f) n=100 replay: mean weight-layer SE / empirical SD ∈ [0.90, 1.10] at every rho; rho=0 coverage ∈ [0.92, 0.97] → ___ [PENDING]
- [ ] (g) kink count (|lambda t_j − 1| < 1e-6) reported; ≤ 1% of fold×source cells → ___ [PENDING]

### 5. Validation results (filled after running)
PENDING — Slurm job 18529696 (`diagnosis/weight_layer/run_weight_layer.sh`), submitted 2026-09-07.

### 6. Decision + rationale
PENDING

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
PENDING — Slurm array 18529697 (`diagnosis/dgp_common_basis/run_dgp_common_basis.sh`, 13 cells), submitted 2026-09-07.

### 6. Decision + rationale
PENDING

# Baseline Formula-to-Code Mapping (FACE-HD)

This note documents how each baseline estimator is implemented and where the key formulas appear in code.

## 1) Target-only AIPW

- Estimand: $\mu_t^a = E_t[Y(a)]$.
- Pseudo-outcome:
  - $\phi_i = \hat m_a(X_i) + I(A_i=a)\{Y_i-\hat m_a(X_i)\}/\hat p_a(X_i)$.
- IF variance:
  - $\hat V_{ot} = n^{-1}\sum_i \hat\varphi_i^2$.
  - $\widehat{\mathrm{Var}}(\hat\mu_t^a)=\hat V_{ot}/n$ (with nuisance-adjusted IF in helper).

Code:
- `estimate_target_only`: [R/estimators_target.R](R/estimators_target.R#L24)
- Cross-fit target-only: [R/estimators_target.R](R/estimators_target.R#L59)
- AIPW IF helper: [R/estimators_helpers.R](R/estimators_helpers.R#L77)

## 2) Sample-size weighted (SS)

- Point estimate:
  - $\hat\mu_{SS}=\sum_j (n_j/N)\hat\mu_j$.
- Heterogeneity-adjusted variance (RE):
  - $\widehat{\mathrm{Var}}_{RE}=\sum_j (n_j/N)^2\{\widehat{\mathrm{Var}}_j+\hat\tau^2\}$.

Code:
- Estimator: [R/comparison_methods.R](R/comparison_methods.R#L37)
- DL heterogeneity (`Q`, $\tau^2$, $I^2$): [R/estimators_helpers.R](R/estimators_helpers.R#L423)

## 3) Inverse-variance weighted (IVW)

- Point estimate:
  - FE: weights $\propto 1/\hat V_j$.
  - RE: weights $\propto 1/(\hat V_j+\hat\tau^2)$ when $\hat\tau^2>0$.
- Variance:
  - RE variance plus pooled IF safeguard:
  - final variance is `max(var_random_effects, var_pooled)`.

Code:
- Estimator: [R/comparison_methods.R](R/comparison_methods.R#L149)

## 4) Tilted-AIPW baseline

- Per-site source estimate uses density-ratio weighted AIPW.
- Adds first-order density-ratio score correction (`adj_alpha`) to source IF.
- Final variance decomposition:
  - `target_variance_component + source_variance_component`.

Code:
- Estimator + decomposition: [R/comparison_methods.R](R/comparison_methods.R#L309)
- Decomposition terms: [R/comparison_methods.R](R/comparison_methods.R#L507-L536)

## 5) Federated-DR baseline

- Source sites:
  - DR weights via exponential tilting against target moments.
  - Site AIPW IF + first-order DR-weight correction.
- Aggregation:
  - IVW point estimate from site-level variances.
  - Variance decomposition consistent with shared-target structure:
    - `var_target_component + var_source_component`.

Code:
- Estimator: [R/comparison_methods.R](R/comparison_methods.R#L553)
- Shared-target/source decomposition: [R/comparison_methods.R](R/comparison_methods.R#L651-L675)

## 6) Pooled-DR baseline

- Pools all observations with site-specific DR weights.
- Uses weighted AIPW IF from helper.
- Adds first-order DR-weight correction per source (source + shared-target parts).
- Final variance uses site-stratified within-group centering:
  - $\widehat{\mathrm{Var}} = N^{-2}\sum_g\sum_{i\in g}(\hat\varphi_i-\bar\varphi_g)^2$.

Code:
- Estimator: [R/comparison_methods.R](R/comparison_methods.R#L690)
- Stacked correction + WSS variance: [R/comparison_methods.R](R/comparison_methods.R#L760-L803)

## 7) Shared helper pieces

- DR lambda CV for density-ratio model:
  - [R/estimators_helpers.R](R/estimators_helpers.R#L475)
- DR weights `exp(-eta)` with normalization/clipping:
  - [R/estimators_helpers.R](R/estimators_helpers.R#L513)
- Weighted AIPW now uses nuisance-adjusted IF helper:
  - [R/estimators_helpers.R](R/estimators_helpers.R#L557)

## 8) Notes on strictness

- Implementations are aligned at first-order IF/asymptotic level.
- Finite-sample exactness is not guaranteed for any semiparametric baseline.
- The code avoids heuristic `dr_inflation` and uses explicit IF-based variance paths.

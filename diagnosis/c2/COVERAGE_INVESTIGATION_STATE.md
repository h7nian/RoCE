# C2 Coverage Investigation — Current State

**Question:** Why is CI coverage a little low in some C2 settings (p=10, K=3/4)?

## Established facts (evidence-backed)

1. **Structure is correct.** method-alignment audit + gated tests pass: FACE aggregation-CV,
   SMMAL two-layer cross-fitting, `docs/main.tex`, and the implementation agree. Not a
   structural CV / cross-fitting bug.

2. **Propensity is NOT the cause.** Channel decomposition of the on-disk
   `score_moment_audit/results/*fold_source_scores.csv`:
   - pure π channel `model_delta_true_pi − model_delta_observed_A` = −0.001…+0.006
     vs required δ ≈ 0.015–0.030. Swapping true π does almost nothing.
   - The prior "true π recovers the moment" was an artifact: `model_delta_true_pi` swaps the
     observed indicator AND the residual `(Y−m̂)→(m_true−m̂)` at once; the residual swap is
     near-circular with `required_delta`.

3. **Signature:** `observed_delta ≈ 0` robustly; same weight against the TRUE residual
   (`mdl_obsA = E_s[ŵ(m_true−ψ̂)]`) ≈ `required`; the gap is
   `out_noise = E_s[ŵ(Y−m_true)] ≈ −required`, which should be ≈0 in population.

4. **Cross-fitting is proper.** `process_source_site` (R/cross_fitting_algorithms.R:141)
   trains nuisances on `setdiff(1:n_folds, c(k1,k2))` and evaluates δ on fold k1. So
   `observed_delta≈0` is a genuine out-of-fold property, not leakage at the k1 level.

5. **Mechanism (code + theory).** The calibrated outcome `fit_unified_outcome_cpp`
   (src/outcome_model.cpp:90,108,162) is a weighted logistic of Y on X with weights
   `w=exp(−Zγ̂)` and an intercept. Its intercept FOC is `Σ w_i(Y_i−ψ̂_i)=0 ⟺ E_s[ŵ(Y−ψ̂)]=0`,
   which **is** the correction δ. `main.tex` eq:alpha_calibrated_loss (line 653) has the same
   FOC (φ includes an intercept, line 237), and the paper says it "enforces U_γ=0". So at the
   calibrated solution δ*→0 and the source estimator collapses to the biased plug-in
   `E_t[ψ̂] ≈ μ₁ − required`. Code faithfully implements the paper → construction-level, not a
   code-vs-paper bug.
   - Algebra: under balancing, estimator bias = `out_noise = E_s[ŵ(Y−m_true)]` exactly.

## Open question (the decisive triage)

Is `out_noise = E_s[ŵ(Y−m_true)] ≈ −required` a **finite-sample higher-order bias**
(vanishes as n grows → coverage recovers, no method change needed) or a **structural** effect
(persists → real method/identification issue)?

At n=5000: bias ≈ −0.019, se ≈ 0.010 → |bias|/se ≈ 2 → coverage ~0.6. Consistent with a
slowly-vanishing higher-order term.

## CORRECTION: it is NOT finite-sample. Coverage DEGRADES with n.

The 4-seed n-scaling (job `9926085`) looked like it shrank, but those 4 seeds were
**selected for high bias at n=5000** → regression to the mean. The unbiased 50-seed coverage
study (job `9927019`, 200/200 done) shows the opposite:

| K | n=5000 | n=20000 |
|---|--------|---------|
| 3 | **0.94** (bias −0.0065, se 0.0103) | **0.84** (bias −0.0041, se 0.0052) |
| 4 | **0.94** (bias −0.0080, se 0.0104) | **0.88** (bias −0.0040, se 0.0052) |

bias shrinks ~1.6× while se shrinks ~2× over a 4× n increase → normalized bias persists/grows
→ coverage degrades toward ~0.85. **This is a first-order / structural bias, not benign
finite-sample.**

## Root cause (strong evidence): the estimated density-ratio tilt γ̂

- **Oracle γ fixes it.** `oracle_gamma_p10_probe`: oracle (true γ) DR has bias ≈ −0.004,
  coverage 1.0 at n=5000, vs target-only bias −0.023 / coverage 0.625.
- **`dr_cv_scale_patch_probe` residual_balance:** with TRUE γ (fitted α), source bias goes
  −0.0159 → **+0.0013** (eliminated). Changing γ regularization moves it (×10 λ → −0.0092).
- So under outcome misspecification (C2), consistency rides on the empirical balancing
  identity `E_s[ŵ·g]=E_t[g]`, which fails because γ̂ is mis-estimated.
- **Mechanism in code:** density-ratio CV validation loss (`src/cv_utils.h`
  `density_ratio_val_loss`, line 386) divides the exp-term by `n_val = n_treated_val`
  (CV folds TREATED-only indices: `density_ratio.cpp` `filter_treated`, ~line 316), i.e. a
  **treated-arm average**. But the training/population objective normalizes by `/n` (full
  source/total): `density_ratio.cpp:152`, comment line 34. **Scale mismatch in λ_γ
  selection → mis-regularized γ̂.** Matches the prior chat's flagged candidate.
  NB: a prior `c2_dr_cv_scale_patch` attempt "improved but didn't fully solve" — understand
  why before re-patching.

## Running experiment: oracle-γ coverage vs n (job `9929392_[1-240]`)

- `c2_oracle_gamma_scaling.{sh,cmd}` (reuses unmodified `c2_oracle_gamma_probe.R`):
  60 seeds × cells {(K3,n5000),(K3,n20000),(K4,n5000),(K4,n20000)}. Output:
  `diagnosis/c2/oracle_gamma_scaling/`.
- **Decision rule:**
  - oracle coverage stays ~0.95 at n=20000 (no degradation) → estimand/identification sound;
    the bug is γ̂ ESTIMATION (the CV-scale mismatch). Then: fix `cv_utils.h` density-ratio CV
    validation scale (treated-arm → full-source), re-run coverage to confirm, then patch src.
  - oracle coverage ALSO degrades → deeper issue (outcome calibration / estimand); revisit the
    derivation tension (calibrated-α intercept FOC zeroing δ).

## CONCLUSION: estimand is sound; bug is γ̂ CV-scale (job `9929392`, 240/240)

Oracle (true γ) coverage is at/above nominal at ALL n and does NOT degrade:

| method | K3 n5000 | K3 n20000 | K4 n5000 | K4 n20000 |
|--------|----------|-----------|----------|-----------|
| oracle γ (raw/source-cond) | 0.967 | 1.00 | 0.983 | 1.00 |
| fitted γ (run_crossfit)    | 0.94  | 0.84 | 0.94  | 0.88 |

oracle bias ≈ 0.0000 everywhere. → identification/estimand correct; the entire C2
under-coverage is **γ̂ estimation**. Prior CV-scale patch was tested only at n=5000 (where
coverage ≈0.94), so it could not reveal the fix.

## FIX IMPLEMENTED (src) — pending validation

`src/cv_utils.h::density_ratio_val_loss` now takes `source_scale` and scales the validation
exp-term by `n_treated/n` (arm fraction), so the CV validation loss estimates the full-source
population risk `Ẽ_{s_j}[I(A=1) exp(-φᵀγ) ψ']` instead of a treated-arm average. The 3 callers
(`src/density_ratio.cpp` refined/initial/calibrated CV, lines 368/442/527) pass
`static_cast<double>(n_treated)/n`. Matches main.tex (eq:gamma_init, eq:gamma_calibrated_loss)
and the training normalization (`density_ratio_cd_update` divides by full `n`).

**VALIDATED ✓ (job `9930454`, 200/200).** Patched-γ coverage no longer degrades with n:

| K | n | baseline | fixed |
|---|---|----------|-------|
| 3 | 5000 | 0.94 | 0.94 |
| 3 | 20000 | 0.84 | **0.94** |
| 4 | 5000 | 0.94 | 0.96 |
| 4 | 20000 | 0.88 | **0.96** |

n=20000 bias roughly halved (−0.0041→−0.0024 K3; −0.0040→−0.0017 K4). Fix confirmed.

Remaining before deploy/commit: (a) regression check on non-C2 configs (fix changes λ_γ
globally) — paired baseline-vs-fixed coverage for C1/C3/(C4) and C2 p=50; (b) testthat suite
against the fixed build; (c) reinstall to `~/Rlibs` + commit only if (a),(b) clean.

---
(historical) patched RoCE built into isolated lib
`diagnosis/c2/fix_validation/Rlib` (build job `9930453`); fixed-coverage array
`9930454_[1-200]` (50 seeds × {K3,K4}×{n5000,n20000}, same seeds as baseline `9927019`),
output `diagnosis/c2/coverage_scaling_fixed/`.
- PASS if coverage rises to ~0.93–0.95 at BOTH n (esp. n=20000: 0.84→~0.95).
- Then: re-run a non-C2 config (e.g. C1/C3) to check no regression, run testthat suite,
  and only then reinstall to `~/Rlibs` / commit.
- If only partial: also align the CV training normalization (train_scale) and re-test.

## Regression + honest caveats (in progress)

Paired same-seed (1000-base) C2, baseline→fixed: bias K3 n20000 −0.0041→−0.0024,
K4 n20000 −0.0040→−0.0017 (~half); sd(bias)/se 1.14→1.05, 1.06→1.00; coverage
0.84→0.94, 0.88→0.96. Fix is real and effective for C2.

CAVEAT (do not over-claim): config_coverage_baseline (3000-base, 40 seeds) shows ALL
configs incl. C1 (bias≈0) at coverage ~0.875–0.925, with sd(bias)/mean(se) ≈ 1.07–1.16 at
n=20000 → a **mild ~5–10% SE under-estimate at large n, general across configs** (separate
from the γ bias). BUT at 40–50 seeds sd/se carries ±0.11 noise, so this is borderline — needs
more seeds to confirm it is real before any variance-formula change. Coverage absolute levels
are seed-noisy (±0.05); the fix EFFECT is established via the paired same-seed comparison.

Pending: config_coverage_fixed (3000-base, all C1–C4) for the paired regression (confirm
C1/C3 do not regress, C2 improves); baseline testthat `9936638` to confirm the 2 smoke-test
failures (degenerate folds in target-only/comparison paths, NOT density-ratio) are pre-existing.

## Constraints
- No R in the login terminal; submit SLURM `.sh` jobs only (public partitions).
- Work in `diagnosis/` + `tests/`; modify `R/`/`src/` only after a probe gives significant
  evidence (now met for the γ̂ CV-scale fix; still validating before deploy/commit).

## 2026-05-29 — Two lever probes (both NEGATIVE): M_tau truncation and lambda-grid selection

### M_tau truncation (job 10001706, 60 tasks, C2 p50 n=5000 K=3, 20 seeds x M_tau in {3,4.4,10})
Probe: diagnosis/c2/c2_score_moment_audit.R reads env C2_M_TAU and passes M_tau to run_crossfit.
Result: the calibrated-loss truncation T(phi^T alpha) is **INERT at C2 p50**. Same seed ->
BYTE-IDENTICAL estimate (e.g. 0.5807274535) for ALL of M_tau in {3,4.4,10} -> |phi^T alpha| < 3
never clips. Aggregate (20 seeds each, all three M_tau identical):
  mean_bias=-0.0075  rmse=0.0119  mean_se=0.0104  emp_sd=0.0093  coverage=0.900  bias/se=0.71.
=> SMMAL-style tighter truncation (2M~=4.4) does NOT reduce the C2 p50 bias or lift coverage.
   M_tau is not the missing lever. (Plumbing verified: 20 logs each at M_tau=3/4.4/10.)

### lambda_gamma grid+selection: RoCE dense lambda.min vs RCAL 0.5^k+tune.cut (job 10018630)
Confirmed from source that RCAL::glm.regu.cv (what SMMAL used) selects lambda by
which.min(out-of-fold "cal" ENTROPY) on a coarse grid lmax*0.5^(0:10) (nrho=11, tune.fac=0.5,
tune.cut=TRUE). RoCE uses the SAME entropy loss + lambda.min, on a dense 100-pt grid down to
1e-4*lmax (n>p). New diagnosis script diagnosis/c2/c2_lambda_grid_compare.R installs a RUNTIME
namespace patch on select_lambda_cv_calibrated_density_ratio_cpp (NO source edits) to capture the
live calibrated-gamma CV curve from run_crossfit, then emulates RCAL's pick on the same curve.
24 seeds x 30 curves = 720 curves, C2 p50 n=5000 kf=10:
  ratio RoCE/RCAL: median 0.991, mean 0.985, [q05 0.69, q95 1.31]
  |log2 ratio|: median 0.21, mean 0.23, MAX 0.745  -> 100% within ONE RCAL 2x grid step
  frac RoCE lambda SMALLER: 0.508 (SYMMETRIC -> NO systematic under-regularization)
  frac RoCE at grid floor: 0.000  -> the 1e-4-vs-9.8e-4 floor difference is moot; optimum interior
  entropy_gap (RCAL-RoCE): median 2.3e-4, max 2.9e-3 -> RCAL coarse grid loses ~nothing
=> RoCE's lambda.min and RCAL's lambda.min are STATISTICALLY EQUIVALENT at C2 p50. RoCE
   faithfully reproduces SMMAL/RCAL not just in formula but in the SELECTED VALUE. The grid-design
   differences (denser grid, lower floor, no tune.cut) do NOT cause under-regularization here.
   Perf note: this CV is slow only at small n (n/p small: cv_secs~299s at n=263) but ~2-3s at
   n=5000 (n/p~25); run smokes at the REAL n, not shrunken n.

### CONSOLIDATED CONCLUSION on C2 p50 ~0.90 coverage
Three independent probes now converge: (1) n-scaling -> bias shrinks with n (->0.95 at n=20000);
(2) M_tau truncation INERT; (3) lambda-selection EQUIVALENT to the SMMAL/RCAL reference and never
at floor. => the C2 p50 under-coverage is the IRREDUCIBLE finite-sample high-dim nuisance-estimation
bias (bias/se~0.7), NOT a fixable truncation or lambda-selection artifact. Any further improvement
requires a DEPARTURE from SMMAL (targeted/MSE-aware lambda, or higher-order debiasing), not a
re-tuning of the existing CV.

## 2026-05-30 — Direction B: targeted-lambda (balancing-MSE) prototype — NEGATIVE (it HURTS)

Built a SMMAL-departing prototype: select lambda_gamma by minimizing the OUT-OF-FOLD BALANCING-MOMENT MSE
||U_gamma||^2, U_gamma = E_t[grad_alpha psi] - E_s[I(A=1) e^{-phi^T gamma} psi' phi] (the first-order bias driver,
= fold_gamma_score), instead of the density-ratio ENTROPY. All pieces adversarially verified (workflows):
  - c2_targeted_lambda_probe.R: runtime namespace patch on select_lambda_cv_calibrated_density_ratio_cpp captures
    the live entropy curve AND computes a warm-started K_bal=3 balancing-MSE CV curve on the same calibration rows.
    GATE (job 10028635, 2 seeds): lambda_balancing DIFFERS materially from lambda_entropy -- median(lam_ent/lam_bal)
    ~0.72-0.83 => balancing wants ~20-40% MORE regularization (NOT symmetric like RoCE-vs-RCAL). Gate passed -> proceed.
  - c2_targeted_lambda_coverage.R: OVERRIDE patch forces the estimator's lambda_gamma = lambda_balancing (sets
    res$lambda_min/lambda_1se/best_lambda + idx), paired against baseline (entropy) on identical data/folds/seed
    (only lambda differs). Verified against R/model_fitting.R:102-126 that overriding lambda_min+lambda_1se is
    sufficient for rule="min". Balancing search on a stride-4 subgrid (always incl. idx_entropy).

RESULT (job 10051686, n=40 matched seeds, C2 p50 n=5000 kf=10):
  coverage   BASE(entropy) 0.925  /  TGT(balancing) 0.875   (Delta -0.050)
  mean bias  -0.00822 / -0.00945 ; mean|bias| 0.01042 / 0.01097 ; rmse 0.01237 / 0.01291
  mean se     0.01049 / 0.01035  (only ~1.3% variance reduction)
  median lam_ent/lam_bal = 0.737 (targeted ~36% MORE regularized)
  paired Delta bias = bias_tgt - bias_base: mean -0.00123, se 0.00024, t = -5.11; targeted MORE biased 28/40.
=> TARGETED-LAMBDA HURTS. The OOF balancing MSE = bias^2 + VARIANCE of the moment; at finite n the variance term
   dominates, so it picks larger lambda; that buys ~1% se but a significant BIAS increase (t=-5.1). Since C2 p50 is
   BIAS-LIMITED, the trade is strictly counterproductive.

CONSOLIDATED (closes the lambda question both directions): entropy-CV lambda is already NEAR BIAS-OPTIMAL for the
estimator. Smaller lambda (earlier dr_lambda_coverage / residual_balance probes) -> coverage flat as variance rises;
larger lambda (this balancing-targeted run) -> bias rises, coverage falls. So C2 p50 under-coverage is NOT improvable
by ANY lambda_gamma re-selection (entropy, RCAL-equivalent, or balancing/MSE-targeted). It is the irreducible
finite-sample high-dim nuisance-estimation bias; only n (or a higher-order/debiased correction, a separate research
direction) can move it. M_tau truncation also inert. Recommendation: keep the current CV-on-entropy lambda.min (faithful
to SMMAL/RCAL and near bias-optimal); do NOT adopt balancing/MSE-targeted lambda.

## 2026-06-01 — RULES UPDATE / CORRECTIONS (supersedes earlier conclusions; read as current truth)

### (A) "kf=10 too fine" / "n ∝ K_f" operating-range rule is from OLD code
- diagnosis/bias/README's mechanism (per-fold calibration too small -> exp(Zγ) blow-up ->
  colMeans PARAMETER-AVERAGING across k2 -> O(1) Jensen gap) was measured on the INITIAL-commit
  code (gamma_final <- colMeans(per-k2 γ), 98040a2).
- CURRENT working tree REFACTORED to FOLD-SUMMED/stacked calibration
  (R/cross_fitting_algorithms.R:195-216: .stack_fold_field + .make_plugin_block_design + ONE
  fit_unified_density_ratio on the stacked (K_f-1) folds). `git log -S "Z_cal_stack"` => uncommitted.
- => per-fold-too-small mechanism largely GONE; calibration uses ~(K_f-1)/K_f of the site even at
  kf=10 (LARGE). "kf=10 too fine" and "n ∝ K_f" are NOT established for current code.
- RE-TESTED NOW: diagnosis/face_probe (10185918 dev=0; 10186618 deviation), kf=5 vs 10 at fixed n.
  cov ~equal => rule obsolete; kf=10 undercover => still holds.
- WHY n∝K_f would hold: binding quantity is per-fold n_calib_arm = n_site_arm/K_f (shrinks with K_f).
  Standard 1-layer CV does NOT (more folds = more training data); it was a quirk of the 2-layer +
  nonlinear-exp-nuisance + parameter-averaging construction, not a general CV fact.

### (B) The aggregation is FACE-faithful and NOT the C2 culprit
- RoCE aggregation = FACE eq 9/11: anchor (target-only) + Σ η(Δ_ks-Δ_T); penalty weight
  (mu_ot-estimate)^2 = (Δ_T-Δ_ks)^2 (weight_optimization.cpp:14; cross_fitting_aggregation.R:410).
  Inner-CV objective is VARIANCE-only (.validation_aggregation_objective, lambda=0).
- DECISIVE: aggregate = affine combo, weights sum to 1 => COMMON-MODE bias b (shared by target
  anchor + sources) is INVARIANT to η. NO weight objective (variance/MSE/penalty) can remove it.
- C2 bias = common-mode high-dim nuisance bias (target_only is ALSO biased, even more than RoCE).
  Coverage gap target_only(~0.94) vs RoCE(~0.90) is PURE bias/SE (RoCE halves SE, bias barely
  shrinks). The idiosyncratic per-source part IS averaged out -> K-effect (cov rises with K).
- => "fix the aggregation objective / add penalty to the CV" does NOT help C2. (C3 shows RoCE robust,
  so variance-CV adequately protects against negative transfer.)

### (C) targeted-λ (balancing-MSE) — CLOSED, it HURTS
- B experiment (10051686, 40 paired seeds, C2 p50): entropy-λ cov 0.925 vs balancing-MSE-λ 0.875,
  paired Δbias t=-5.11. More regularization -> more bias. Entropy-CV (=SMMAL/RCAL) is near
  bias-optimal; do NOT change λ selection.

### (D) C2 vs C3/C4 — where RoCE shines (the selling point)
- C2 (DR-consistent sources): mild common-mode bias, bias-limited undercoverage; only n / kf /
  higher-order debiasing can move it.
- C3/C4 (PS/site misspecified -> INCONSISTENT sources): source bias is ASYMPTOTIC (does NOT vanish);
  as n grows baselines (federated_dr, tilted_aipw) CRASH (cov 0.77->0.62 at n 5k->10k), RoCE stays
  robust (~0.93) by shrinking divergent sources out (adaptive penalty). RoCE's headline result.

### (E) DGP rules (for the paper)
- roce (current core DGP) != FACE DGP. roce = binary, model-misspec C1/C2/C3, n split across sites,
  NO ATE deviation -> NEVER tests negative transfer (RoCE's headline).
- FACE DGP = continuous (we use binary on purpose, matching the application), quadratic-truth-fit-linear,
  skewed-normal covariates, 8-level ATE deviation, K=9 n_k=200. RoCE's `face` DGP already SPARSE
  (min(4,p) active, constants.R:338) -> high-dim = just add noise covariates (no code change).
- MAIN LINE = FACE deviation DGP, high-dim sparse, binary, in-regime. n_k MODEST (1000-2000,
  federated-realistic) NOT huge; n=20000 was TOO LARGE (low-p wasteful + undermines federated motivation;
  high-p still out of regime). kf MATCHED to dimension (being tested). per-site n_k = n_total/(K+1).

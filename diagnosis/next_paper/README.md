# Next-paper research workspace

This directory develops weak-deviation and growing-source-count inference.
The current paper's estimator and frozen Slurm campaigns remain separate.
Numerical references, fitted-score diagnostics and conditional derivations
must not be described as a completed full-method theorem.

## Joint contrast and conditional exact calibration (3 October 2026, v18)

NEXT-PAPER JOINT CONTRAST AND EXACT CALIBRATION (2026-10-03, v18): completed
19,000 distinct setting-specific repetitions at 1000 patients/site. Paired
reanalyses are not recounted. Joint noise/nuisance protection improves the
original fitted-score method, but it remains limited against a stronger
unsplit, range-adaptive target-only reference.

A separate four-stratum conditional calibration benchmark cancels the shared
outcome prediction exactly. On 500 fresh weak-bias repetitions, K128/512/2048
protected lengths are .1923/.1646/.1493 versus .3083 for equally protected
unsmoothed target-only. Conditional-outcome and target-composition uncertainty
are both protected. It requires cellwise source validity and sufficient cell
counts; it is not the general high-dimensional production estimator.

On 1000 fresh K512 strong-bias repetitions, projecting the target estimate
onto the conditional interval reduces RMSE .0439 to .0249 (target .0312),
but midpoint remains better under weak bias and arm cancellation. Both fixed
point rules are retained; no scenario-specific rule is selected after seeing
truth. On 500 fresh weak-shift repetitions without shared target heterogeneity,
naive-pooling coverage is .382/.146/.004 as K128/512/2048; protected coverage is1.
All-valid naive coverage remains near95%; its estimated variance is independently
audited against the analytic DGP variance. Forty-five focused tests pass.

Published to the dedicated RoCE-K Overleaf project at `b88210b178d50df6c972c1425d8268806b07f615`;
entries `main.tex` and `experiments_20261003b.tex`.
Design v5 adds joint_contrast.tex and conditional_standardization.tex;
legacy_design_v4.tex preserves the prior entry. Artifacts and full qualifications:
implementation/next_paper/v18/joint_contrast_20261003_v1/report/REPORT_zh.txt.
Current-paper estimators and historical results are unchanged.

Entry points in `residual_compatibility/`: `run_joint_contrast.py`,
`run_joint_followup.py`, `run_conditional_standardization.py`,
`audit_conditional_balance.py`, and `review_joint_contrast.py`. The exact-balance
runner accepts `--fresh`, `--point-check`, or `--weak-stress`; these select
prespecified independent seed panels. `--repeats` and `--workers` control scale.
`source_cohort="heldout"|"all"`, `dispersion_mode="arms"|"contrast"|"adaptive"`,
and `point_rule="midpoint"|"target_projection"` expose the inference comparisons.
Default midpoint behavior is retained. Adaptive dispersion splits its error
budget before comparing the two certificates. Cells lacking variance support
fall back to protected target-only; repetitions are never silently dropped.

The 19,000 total comprises 4500 primary, 6000 independent joint-score checks,
2500 fresh conditional-calibration, 3000 point-rule, and 3000 weak-stress draws.
Oracle truth/variance are used only for evaluation/audits, not the feasible
exact-calibration intervals. All saved reported population intervals target .10;
the finite-target conditional truth is an intermediate diagnostic only.

## Conditional remainder control and design reuse (3 October 2026, v17)

Fetched the new RoCE-K project at `1e15b0e`. The research design is now version
4, adding `roce_k/sections/training_average.tex`; version 3 remains available
as `legacy_hybrid_v3.tex`. Current-paper RoCE/ENAR defaults are unchanged.

Completed 9,700 distinct setting-specific repetitions: 16 x300 main settings,
4 x1000 independent confirmations, and 3 x300 exploratory K2048 settings.
The source-design reuse comparison uses the same 4,800 main datasets and is
not counted again. Every site has 1,000 patients. Thirty-three focused tests
pass; representative final-source reproductions match saved results exactly.

The conditional-outcome bound controls the complete drift in a fixed but
unknown valid subset. It requires saturated weights fitted from X,A without
the same training outcomes. The empirical-design certificate needs no supplied
source probability minimum. A separate binomial MGF option assumes a declared
joint cell/arm probability minimum; .02/.04/.06 are prespecified sensitivities.
The direct-dispersion mode can reuse held-out designs in its inference sample,
so sources use 500 training +500 evaluation rather than 500+250+250. The
existing Fourier implementation rejects this reuse option.

Independent-confirmation K512 weak-bias lengths are .5253 (conditional design)
and .4727 (MGF with .04 minimum), versus about .697 for the equally improved
two-fold target AIPW certificate. However, an unsplit saturated target-only
finite-sample reference has length .3951 and is still shorter. At K2048 weak
bias the lengths are .4422 and .3927, respectively; only the supplementary MGF
option approaches the stratified target reference. Observed coverage is 1
throughout. Strong-bias midpoint RMSE remains worse than target-only; coverage
must not be confused with robust point estimation. No general high-dimensional
borrowing advantage is established.

Artifacts: `/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v17/training_average_20261003_v1/`.
Published to the dedicated RoCE-K project at `695fb9b1c636ee55e125e8c18c75b099b908520f`;
entries `main.tex` and `experiments_20261003.tex`.

Read `report/REPORT_zh.txt`, `compiled/training_note.pdf`, and `compiled/main.pdf`.
Source entry points in `residual_compatibility/` are `run_training_average.py`,
`audit_stratified_reference.py`, and `review_training_average.py`. Controls:
`--reuse-source-design`, `--confirmation`, `--large-k`, `--repeats`, `--workers`.
The conditional certificate selects `source_calibration="empirical_design"`
or `"training_mgf"`; saved method ID `conditional_pilot` is retained for
reproducibility and denotes the empirical-design method, including reused designs.

## Aggregate calibration and component ablation (2 October 2026, v16)

The research design is now version 3 in `roce_k/sections/hybrid.tex`. It defines
Fourier/dispersion source constraints, the three-dimensional target region,
error allocations, shorter fallback, and coverage reduction in one place.
The earlier moment design is available through `roce_k/legacy_moments.tex`;
its Gaussian comparison now distinguishes known-variance `D-circ` from
random-variance-corrected `D-hat`. The joint GLS legacy entry is preserved.

Completed 18 binary settings x 1,000 and 12 fresh finite-stratum fitted settings
x 300 = 21,600 setting-specific repetitions. Direct patient-level exponential
dispersion and pooled empirical Bernstein avoid individual source variance
certificates in the dispersion module. A pilot calibration-moment certificate
protects signed means and squared nuisance drifts on unknown valid subsets.
The fitted model is still a four-stratum diagnostic, not high-dimensional RoCE.

At K=2,048, weak bias, and target heterogeneity .3, the binary interval-length
ratio to ordinary target-normal falls from 1.736 (previous hybrid) to 1.079
(direct hybrid) or 1.066 (direct dispersion). Without target heterogeneity,
direct dispersion reaches .301. All observed coverages remain 1 in this study.
For fitted K=128, weak bias and heterogeneity .3, source budgets fall from
.285/.263 to .111/.097 and length from 1.110 to .872. However, the equally
improved full-target certificate has length .831. No net fitted borrowing
advantage is established. Fourier does not add a net length benefit under
these protected bounded-score calibrations; retain it as an explicit ablation.

Artifacts: `/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v16/unified_inference_20261002_v1/`.
Published to the dedicated RoCE-K Overleaf project at `1e15b0ec564b8ad94aa0c762f30faea53d578838`.
Main entry: `main.tex`; paired experiment report: `experiments_20261002b.tex`.
The current RoCE/ENAR project was not modified.

Read `report/REPORT_zh.txt`, `compiled/concentration_note.pdf`, and the updated
`compiled/main.pdf`. `paired_target_audit/` reuses all 3,600 fitted samples;
it does not add Monte Carlo repetitions. Its certificate audit supersedes the
raw strict-zero diagnostic, which flagged roundoff in an exactly zero drift.
No intervals or coverage estimates were changed by that diagnostic correction.

Entry points: `run_concentrated.py --profile binary|fitted --repeats N --workers W`,
`audit_concentrated.py`, and `review_concentrated.py` under
`residual_compatibility/`. `infer_repair` accepts `source_mode` = `hybrid`,
`fourier`, or `dispersion`, and `dispersion_calibration` = `cantelli` or
`bounded_exponential`. Defaults preserve historical research behavior.
The new `source_squared_bias_budget` accompanies the signed subset-mean budget;
without it, `source_bias_budget` retains its uniform-coordinate meaning.
Production RoCE/ENAR defaults and earlier results are unchanged.

## Latest follow-up (2 October 2026)

The fallback/calibration follow-up is complete: 18 Gaussian settings with
2,000 repetitions, 18 binary-score settings with 1,000, and 12 fresh
saturated-nuisance settings with 300. The 57,600 setting-specific repetitions
all respect 1,000 patients per site. A paired nuisance-budget refinement uses
the same observations and is not counted as new independent data. Ninety
archived actual RoCE fits and 40,000 independent target-reference repetitions
were audited separately.

The fixed-budget shorter fallback has a checked pointwise source-interval
length envelope and substantially reduces the earlier fallback-driven RMSE.
Decreasing the target error allowance can widen intervals when shared target
noise is present. Bernstein target regions improve the binary-score result,
but estimated-variance certificates and feasible nuisance budgets remain costly.
Paying the shared target nuisance event once removes an unnecessary K cost;
the fully feasible fitted interval is still about 1.10--1.11 times a target-only
interval with comparable finite-sample protection in the illustrated setting.
It is not yet an efficient full fitted procedure.

New components are `repairs.py`, `run_repairs.py`, `fitted_strata.py`,
`run_fitted_repairs.py`, `export_fitted_packets.R`, `test_repairs.py`,
`audit_target_reference.py`, `review_repairs.py` and `repair_note.tex` under
`residual_compatibility/`. The saturated fitting model is a low-dimensional
diagnostic, not the production sparse/fold-summed RoCE program. The archived
RoCE cases are single-fit diagnostics, not new coverage repetitions.

Artifacts are under
`/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v15/calibration_fallback_20261002_v1/`.
Read `report/REPORT_zh.txt` or `compiled/repair_note.pdf`. The final fitted
panel is `fitted_shared_budget/`; the first paired version remains in `fitted/`.
The standalone report is also published in the RoCE-K Overleaf project as
`experiments_20261002.tex`, verified commit
`c2c61db0b7d453616e2a0231ea4095acf347d6a3`. Existing manuscript sources are unchanged.

## Latest experiments (1 October 2026)

`residual_compatibility/` implements the subsequent three-component target
decomposition `(prediction contrast, treated residual, control residual)`.
The research prototype compares Fourier source constraints, dispersion,
coordinate-count search, pooling and target-only references. Fourier sets
are inverted analytically and all residual-set branches are retained during
target-ellipsoid projection. It is not a Guo/RIFL implementation or a fitted
high-dimensional RoCE release.

The reported study has 114 settings: 90 settings with 1,000 repetitions,
18 independent-seed confirmation settings with 2,000 repetitions, and six
actual-minority-valid boundary settings with 2,000 repetitions. Source
counts range from 8 to 2,048. Every reported site has 1,000 patients; the
variance-pilot panel uses 250 pilot plus 750 evaluation patients. Its original
six extra-pilot runs are retained as auxiliary diagnostics and superseded in
the final reported panel.

Known-Gaussian results favor dispersion for weak shifts and Fourier for
strong shifts. A prespecified split of the source error budget lets their
intersection combine these advantages in independent confirmation. Coverage
is conservative. Estimated-variance and bounded-score protection have visible
length costs. Rare target-interval fallbacks dominate some midpoint RMSEs,
and a fixed fallback probability can obstruct expected-length improvement
in K; this needs a separate repair before a full efficiency theorem.

Entry points are `experiments.py`, `test_inference.py`, `review.py` and
`experiment_note.tex` in that directory. Nine tests pass; the final audit
recomputes all 984 reported method-setting rows from saved repeat records.
Results and immutable code snapshots are under
`/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v14/residual_compatibility_20261001_v1/`.
Start with `report_final/REPORT_zh.txt` or `compiled/experiment_note.pdf`.
The uncorrected preliminary report is retained in `report/`; use
`report_final/` for the final patient-budget-corrected results.

## Previous method design (1 October 2026)

The version 2 design was in `roce_k/main.tex`, with method and
conditional proofs in `roce_k/sections/compatibility.tex` and verified references in
`roce_k/references.bib`. Its dedicated Overleaf project is
https://www.overleaf.com/project/6abed6440abbc7c132b9611d.
It is separate from the current RoCE/ENAR manuscript.

The design prioritizes a confidence region from common-outcome source means,
noise-corrected source dispersion, arm-specific valid-count assumptions and
a target anchor. The final joint mean region has at most four coordinates
and projects onto TATE. It is a conservative moment relaxation of count
compatibility, not Guo's full searching-and-sampling procedure. One-dimensional
borrowing is an optional point-estimation family. The previous joint GLS
design is preserved as `roce_k/legacy_gls.tex` with `sections/design.tex`.

An independent pilot supplies variance upper bounds; simultaneous valid-score
bias budgets remain inputs requiring justification. The first proof retains
independent training, pilot and evaluation roles. Uniform coverage for growing
K, useful constants and efficient full cross-fitting need their stated
conditions; fixed-dimensional final optimization does not establish them.

`check_center_dispersion.R` passes 990 assertions, including 4,096 exact
non-Gaussian score configurations and 100 joint-noise-region witnesses.
The ideal Gaussian radius diagnostic confirms that Cantelli dispersion
protection can overwhelm borrowing gains even without nuisance uncertainty.
These are formula and algebra diagnostics, not complete-method coverage or
expected-length experiments. Tighter dispersion calibration and practical
bias budgets take priority over larger fitted simulation campaigns.

`check_joint_quadratic_bound.R` passes 696 deterministic assertions,
including 60 matrix cases and all 262,144 configurations of a small binary
observed-data example. These are conditional mathematical checks, not fitted
coverage experiments. The next priority is useful nuisance-bias and
covariance bounds, followed by a complete fitted-method validation.

Current artifacts, source snapshots, checks and PDF builds are under
`/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v13/method_design_overleaf_20261001_v2/`.
The version-1 manuscript and its 696 checks remain under
`implementation/next_paper/v12/method_design_overleaf_20261001_v1/`.
Earlier v1--v12 results and designs are retained as research history. The
older Gaussian constructions below are references and diagnostics, not the
new design's inference justification.

## Code and checks

| Component | Implementation | Verification |
|---|---|---|
| Known-Gaussian confidence-set reference | `source_confidence_sets.R` | `test_source_confidence_sets.R` |
| Gaussian variance-region calibration | `variance_region_calibration.R` | `test_variance_region_calibration.R` |
| Actual candidate means and score covariances | `candidate_score_summary.R` | `test_candidate_score_summary.R` |
| Joint arm covariance and pairwise TATE reference | `arm_candidate_score_summary.R`, `arm_pair_confidence_sets.R` | `test_arm_pair_confidence.R` |
| Bias-dispersion geometry | `dispersion_bias_bound.md` | `check_dispersion_geometry.R` |
| Gaussian dispersion intervals | `dispersion_bias_intervals.R` | `test_dispersion_bias_intervals.R` |
| Biased-source population orthogonality | `probe_biased_source_orthogonality.R` | Quadrature and finite-difference checks |
| Target-local IPW projection | `target_common_projection.R` | `test_target_common_projection.R` |
| Gaussian experiment settings | `gaussian_validation_settings.R` | `test_gaussian_validation_settings.R` |

Source `source_confidence_sets.R` before the variance-region module. Source
`candidate_score_summary.R` before the cross-fitted target-projection wrapper.
The experiment and diagnostic drivers load these dependencies explicitly.
All drivers require output paths under `/scratch.global/zhan9381/FACE-HD/`.
Frozen copies, code hashes, logs and results are versioned under scratch
`implementation/next_paper/v1` through `v15`.

`diagnose_candidate_scores.R` reads a fixed first-seed plan and reconstructs
candidate estimates from completed artifacts. It checks data hashes, fold
coverage, target variance and source-message identities. Missing completed
fits stay visible in its availability report. These are structural checks,
not an MC coverage estimate. The comparison-source labels come from the DGP
for diagnosis; they are not supplied to a proposed deployment algorithm.

`diagnose_target_projection.R` supports `main` and `source_count` profiles.
It reuses held-out-fold propensity fits and compares IPW with ordinary outcome
projections. Ordinary fitting is a diagnostic comparator: under outcome
misspecification its population projection need not match the source limit.

## Inference and communication scope

The extracted covariance is an empirical score covariance. It does not add an
unproved nuisance-remainder correction. The default output retains moments;
individual diagnostic score records require `keep_scores=TRUE`.

Source residual means and variances are computed over each entire site fold,
including the treatment indicators already in the residual scores. They must
not be normalized by the arm-specific sample count. Pooling includes sample
counts and between-fold mean variation. The next scalar-TATE interface needs
those source moments and the outcome coefficients used at the target; the
current arm-specific algorithm also retains the cross-arm covariance.

## Theory and remaining work

Start with `theory.md`, `shared_target_structure.md`,
`variance_uncertainty.md` and `uniform_calibrated_score_remainder.md`.
`half_valid_lower_bound.md` and `majority_oracle_adaptation_bound.md` concern
uniform local efficiency in Gaussian models. Their observed-data embedding
and sharpness remain open. `literature_scope_audit.md` records the comparison
with Guo's JRSSB papers and published RIFL; no novelty claim is made.

Next tasks are fitted-score variance inference, uniform nuisance rates,
arm-specific validity and actual comparison with RIFL. A candidate-only
fitting path could omit computations used solely to learn aggregation weights,
but must first match saved outer-fold fits exactly. The current paper uses
two-layer fitting; its three-layer version remains a diagnostic comparator.
`arm_pair_confidence_reference.md` describes a next reference for different
valid-source sets in the two arms. Its pair dependence and marginal variance
requirements are explicit. The implementation and107 saved-fit moment
identities pass; fitted-method inference validation remains open.
The first24x1000 known-Gaussian study is conservative and generally produces
longer intervals than ordinary target-only inference. It is a validity
reference, not a completed efficiency solution. See
`arm_pair_inference_reduction.md` for the conditional growing-K error bound.

The next candidate is a bias allowance based on source dispersion after
accounting for Gaussian sampling noise. `dispersion_bias_bound.md` records
the derived geometry, its established Q-statistic connection, and material
limitations: Gaussian approximation is needed for biased sources too;
fixed heterogeneity can keep its interval wide; growing-K quadratic
statistics may require stronger nuisance rates. The geometry passed280
random subset checks. The Gaussian interval is implemented and its59 tests
pass. A24x1000 paired comparison still finds it wider than a baseline that
uses the same guaranteed-valid control information. The fitted-method theorem
and a general efficiency solution remain open.
`biased_source_initialization_audit.md` documents the separate issue of
regularity at biased source means under one-round target OR initialization.

# Next-paper research workspace

This directory develops weak-deviation and growing-source-count inference.
The current paper's estimator and frozen Slurm campaigns remain separate.
Numerical references, fitted-score diagnostics and conditional derivations
must not be described as a completed full-method theorem.

## Local precision and safe design-branch budgets (4 October 2026, v22)

NEXT-PAPER LOCAL WEAK BIAS AND DESIGN BUDGETS (2026-10-04, v22): completed
4,800 fresh local-model setting-repetitions and paired reanalysis of 2,800 v21
repetitions. These are not 7,800 independent new datasets. The linear model is
the primary paper branch; arbitrary invalid means remain a broader extension.

Derived an exact Bernoulli all-valid expected-length lower bound for intervals
honest over a local weak-bias class. Count-conditioned priors enforce the valid
source guarantee; target information is included. The construction covers
both fixed treatment-arm samples and an i.i.d. randomized population-TATE
submodel. For comparable sizes it gives the necessary n^(-1/2) K^(-1/4) scale.
Matching the full population upper rate additionally requires constant CATE,
fixed dimension, stable designs, K=o(n^2), and negligible bad-design cost.
This is not optimality for arbitrary nonlinear invalid means or all scenarios.

A common coefficient/composition event now permits safe conditional budget
reallocation on design-only abstention. With shared composition .0225 and
coefficient .005, the remaining .0225 is conditional: use all of it for target
noise on abstention, or split it among source/target/dispersion when borrowing.
Abstention intervals exactly equal full-budget protected target-only. Legacy
budgets and conditional-only reallocation remain explicit controls.

Paired p4,K2048,weak-bias population length improves .112275 -> .104327;
ordinary target Wald is .122697. Observed protected coverage remains1.000.
The p10 and heterogeneous-CATE controls still exceed Wald length. In the
n_s120,p45,K64 poor-overlap setting, design-only precision gating returns
protected target exactly (.354644 versus legacy .387492). Weight thresholds
2/4/8 are prespecified sensitivities, not selected from outcome performance.

The fresh exact-intercept panel uses n1000/site and K16--16384. At K16384,
mean-matched weak bias c=.25 gives protected length .0120 and coverage1.000;
naive source Wald coverage is .020 (same-direction: .010). This model declares
constant arm means, so no population-composition enlargement is needed.

70 tests pass. All 2,800 legacy repetitions replay exactly; every applicable
shared-budget design fallback equals the full target interval. Independently
recomputed140 budget rows,294 paired comparisons and96 local rows. Previous
results, v8 manuscript, and production RoCE/ENAR remain unchanged.
Artifacts: implementation/next_paper/v22/weak_bias_budget_20261004_v1/.
Published RoCE-K version 9: 69ea1b0ca4620cd4fc886bfae3280a69b764181e; main.tex and experiments_20261004d.tex.

`design_budget.py` implements the common-event ledger, design-only precision
gate and population wrapper. `budget_policy` supports `legacy`,
`conditional_reallocation`, and `shared_composition`. The selection module adds
`weight_inflation_limit`, preserving the original valid-count bookkeeping.
The precision gate is a design proxy, not an optimality guarantee.

`weak_bias_precision.py` evaluates exact lower bounds, with `experiment` set to
`fixed_arm` or `iid_randomized`; sample-size units are recorded in the output.
`run_weak_bias_precision.py` generates the exact-intercept panel.
`run_budget_followup.py --previous <v21-root>` runs paired confirmation/removal
reanalyses; controls include `--profile`, `--model-mode`, `--repeats`,
`--workers`, and `--run-name`. `review_weak_bias_followup.py` audits endpoints
and generates the tables and figure on scratch. `code_executed/` preserves
byte-identical running sources, while `code_final/` records the verified final
implementation. Main text adds `sections/local_precision.tex` and
`sections/design_budget.tex`; the report is `experiments_20261004d.tex`.

## Main theorem, evaluated removal and precision boundaries (4 October 2026, v21)

NEXT-PAPER MAIN GUARANTEE AND DESIGN REMOVAL (2026-10-04, v21): completed
4,400 outcome setting-repetitions (22 settings x 200) plus exact finite-model
calculations. The dedicated RoCE-K main text now has one algorithm/theorem,
with separate linear, declared-approximation and bounded-invalid guarantees.

Independent-seed confirmation at p=4, K=2048, h=0, weak bias gives length
0.112275 versus normal target 0.122697: paired difference -0.010421,
MCSE 0.000607 (8.49% shorter). Observed coverage is 1.000 versus 0.930;
200 repeats do not establish undercoverage of the normal comparator. The p=10
control remains longer than normal. All four confirmation settings improve
on the equally protected target-only comparator. No oracle range is supplied.

Design-only removal is now an evaluated interval, controlled by design_policy.
Every source quantity is recomputed at K'=K-m and g'_a=max(0,g_a-m), including
unknown-invalid projection budgets. Source n=120, p=45, poor overlap can retain
eligibility without shortening intervals. With n=1000/site, p=4, K=256 and a
prespecified 1/8 rank-deficient valid-source subgroup, removal reduces length
0.3092 to 0.1617 (bounded-invalid: 0.1756; protected target: 0.2825), with
conditional coverage 0.995 and population coverage 1.000.

Precision analysis separates an explicit certificate floor from a new
worst-case conditional-design lower bound. Arbitrary bounded invalid means
permit overlapping source-data mixtures despite different target effects;
the valid-count conditioning penalty and trusted target information are
included. This does not prove all-valid adaptive optimality, a weak-shift-only
lower bound, or an i.i.d. random-design population lower bound. On three-point
support, a prespecified saturated basis removes the variance-centering cost:
at K=262144, eta=.4, conditional length is 0.00547 versus 0.02575 for the broad
linear-basis certificate. The extra finite-support structure is explicit.

62 tests pass; seven complete repetitions replay exactly, and four saved v20
broader-model repetitions retain identical intervals. All 246 method-summary
rows and 72 paired comparisons are recomputed from endpoints. Old results,
v7 manuscript and production RoCE/ENAR remain unchanged.
Artifacts: implementation/next_paper/v21/main_theorem_removal_20261004_v1/.
Published RoCE-K version 8: fe6411e1852fb7d4920570e5071780a4aa27113b; main.tex and experiments_20261004c.tex.

Entry points in `residual_compatibility/`:

- `design_selection.py`: outcome-independent design preparation, count adjustment,
  retained-source evaluation and complete inference using the common engine.
- `run_design_followup.py`: `--profile confirmation|removal`, `--design-policy
  all_required|remove|paired`, `--model-mode linear|bounded_invalid|paired`,
  `--repeats`, `--workers`, and a new `--run-name` on scratch.
- `run_model_boundary.py`: paired finite-support intervals and exact contamination
  overlap calculations; `--repeats` controls the outcome panel.
- `review_design_followup.py`: paired MCSE, exact binomial coverage intervals,
  independent endpoint audits and generated manuscript assets. Pass the new root
  and `--previous` with the v20 root. Report assets stay on scratch/Overleaf.

`sections/main_guarantee.tex` and `sections/precision_boundary.tex` contain the
new theorem and its scope. `legacy_design_v7.tex` preserves the previous layout.
The new report is `experiments_20261004c.tex`; figures/tables are generated in
`report/` and copied to the Overleaf `experiments/20261004c/` directory.

## Target protection, approximation and design growth (4 October 2026,v20)

NEXT-PAPER TARGET PROTECTION AND APPROXIMATION (2026-10-04,v20): completed
7200 outcome setting-specific repetitions (36x200) and7800 design-only draws
(39x200). A4800-repeat paired model-class extension reuses existing draws and
is not counted again. All outcome panels use1000 patients/site.

One observed-target coefficient ellipsoid protects CATE range and empirical
variance; empirical Bernstein reduces composition cost while applying the
same improvement to target-only. At p4,K256,h=.15,length is .1584 versus .2615
for matched universal protection and .2853 for improved protected target-only.
Ordinary target-normal is .1214 there. At p4,K2048,h=0,new length .1119 is
below normal .1228; at p10 it is still longer (.1443 vs .1234). No true range
or constant-effect information is passed to inference.

Model approximation now corrects mean drift,variance-centering bias,and the
quadratic MGF together. A universal bounded-invalid option needs no supplied
nonlinearity amplitude when target/valid means are linear. Its wider class
can cost length and does not automatically retain the linear fourth-root rate.
A fixed-design nonlinear-invalid example at K262144 has old conditional
coverage0 but population coverage1; repaired conditional coverage is1 under
the declared envelope and .990 under bounded-invalid protection. The latter
source band itself covers all repeats; two losses arise in target combination.

The random-design eta=.2 attempt was rejected for negative probabilities at
support corners. The final grid ends at .15,without clipping; old attempts and
verified reusable cells remain archived. Conditional and population intervals,
all radii and design failures are saved. Matrix-Chernoff lower/upper bounds,
leverage control,and expected-width/fallback conditions are stated explicitly.
Design-only source removal is an eligibility diagnostic,not yet a fitted CI.

56 tests pass. All old intervals in24 paired extension cells match exactly.
Artifacts: implementation/next_paper/v20/target_protection_20261004_v1/.
Published RoCE-K version7: 6855f9cdef4d70d49b4a9a1ce326b341c867dd40; main.tex and experiments_20261004b.tex.
Main text follows conditional balance,leave-out correction,dispersion protection,
and target composition; older approaches are retained in the appendix and v6
is preserved separately. Production RoCE/ENAR and old results are unchanged.

New entry points in `residual_compatibility/` are `run_target_protection.py`,
`run_projection_stress.py`, `run_design_stability.py`, and
`review_target_protection.py`. Their controls include `--profile`, `--run-name`,
`--repeats`, `--workers`, and `--bounded-invalid` where applicable.
`population_protection.py` holds the observed-data coefficient/range/variance
certificate and projection-error budgets; `conditional_linear_bands()` remains
the shared conditional inference engine. Budget keys use `mu1`, `mu0`, and
`contrast`; mean-bias arrays order treated before control.

All new target-protection comparisons share the exact conditional bands and
use matched error allocations. The same target certificate improves the
standalone comparator. Universal-range legacy behavior is preserved and checked
against v19. The unknown-amplitude bounded-invalid option requires only bounded
invalid means, but retains linear target/valid means and declared valid counts.
General nonlinear or sparse target fitting remains unfinished.

Read `report/REPORT_zh.txt`, `compiled/main.pdf`, and
`compiled/experiments_20261004b.pdf`. The corrected random-design panel is
`misspecification_v2`; `misspecification_bounded` and `projection_stress_bounded`
are paired model-class extensions. The initial invalid grid is preserved and
excluded from final admissible-setting summaries. No patient draws were dropped
from the final valid settings. Design-only count-adjusted eligibility does not
claim an evaluated precision improvement.

## Continuous linear-balance benchmark (4 October 2026, v19)

NEXT-PAPER CONTINUOUS LINEAR BALANCE (2026-10-04, v19): completed6600
fresh-design repetitions,33 settings x200,n1000/site. The source module now has
a separate continuous-covariate linear-model benchmark: exact signed balance,
unbiased leave-out variance correction,and a patient-level conditional MGF
bound for dispersion. All sources,including invalid ones,must have conditional
means in the fitted linear span. This is not nonlinear or sparse p>n RoCE.

At p20,K256,weak bias,universal-range interval length is .2623 versus .3792
for equally protected target-only; ordinary target-normal is still .1234.
At p200,K64,the corresponding lengths are .3485/.4081/.1535. In-model protected
coverage ranges .995--1. K2048 weak-bias naive coverage is .360,while its
all-valid control is .975 and protected coverage1. Strong-bias projection
improves midpoint RMSE,but midpoint is better for cancellation/weak shifts.
Declared CATE range .5 is a separate stronger-class sensitivity,not a fitted
range. Nonlinear misspecification remains outside the theorem.

51 research tests passed;18 saved repetitions replayed exactly. All231 metric
rows were independently recomputed. Valid-source drift is below1.84e-15,and
mean estimated/oracle source variance ratios range .9990--1.0009. No repeats
were dropped. Old results and current-paper production estimators are preserved.
Artifacts: implementation/next_paper/v19/continuous_balance_20261004_v1/.
Published RoCE-K version6: 0772c6016451401167c8c269158ce4d55fa0470e; entries main.tex and experiments_20261004.tex.
Next: certified target-composition ranges and approximate/nonlinear remainder
control,with the stronger target-only comparator retained.

The new source files in `residual_compatibility/` are `continuous_balance.py`,
`run_continuous_balance.py`, `review_continuous_balance.py`, and their tests.
The runner exposes `--profile main|boundary|stress`, `--repeats`, and `--workers`.
Every source profile conditions on designs,uses all1000 patients,and permits
arbitrarily weak deviations. QR rank/leverage guards are outcome-independent;
unsupported source designs use target-only and unsupported targets use[-1,1].
The main design retains version5 as `legacy_design_v5.tex`.
Read the compiled report or `report/REPORT_zh.txt`; main/normal target comparisons,
paired Monte Carlo errors,and the deliberate nonlinear failure of assumptions
remain explicit. These prototypes are separate from the production R package.

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

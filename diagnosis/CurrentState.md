# FACE-HD current state

Updated 2026-10-03. This file describes the current configuration and evidence.
Earlier state notes are retained under
`/scratch.global/zhan9381/FACE-HD/implementation/r7/state_archives/`.
Live progress belongs in timestamped scratch snapshots and controller records.

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

NEXT-PAPER CONDITIONAL REFINEMENT (2026-10-03, v17): completed9700 distinct
setting-repetitions (16x300 main,4x1000 confirmation,3x300 K2048 exploration).
Reuse comparisons use the same4800 main samples. The fixed-valid-subset scalar
outcome certificate avoids outcome-box corners under design-only weights.
Direct dispersion may reuse design/evaluation observations by a union bound;
Fourier reuse is explicitly rejected. Every site still totals1000 patients.

At K512 weak bias, new length .5253 improves on equally corrected two-fold
target AIPW (~.697), but an unsplit stratified finite-sample target reference
is shorter (.3951). At K2048 weak bias, the design certificate gives .4422;
a supplementary MGF certificate assuming joint probabilities>=.04 gives .3927.
This is conditional finite-stratum progress, not a general high-dimensional
advantage. Strong-bias midpoint bias remains around .05. All observed coverages
are1. Thirty-three focused tests pass; representative output parity is exact.
Artifacts: implementation/next_paper/v17/training_average_20261003_v1/.
Published to the dedicated RoCE-K project at `695fb9b1c636ee55e125e8c18c75b099b908520f`;
entries `main.tex` and `experiments_20261003.tex`.

Read report/REPORT_zh.txt and compiled/training_note.pdf. Main design v4 includes
sections/training_average.tex; legacy_hybrid_v3.tex retains the previous entry.
Current-paper production defaults and historical results remain unchanged.

NEXT-PAPER AGGREGATE CALIBRATION (2026-10-02, v16): completed 18 bounded
oracle-score settings x1000 and 12 fresh four-stratum fitted settings x300,
21600 setting-specific repetitions. Main RoCE-K design is v3 with explicit
Fourier/dispersion/target definition and coverage accounting; v2 remains as
legacy_moments.tex. No production RoCE or ENAR changes.

Direct exponential dispersion avoids individual variance certificates. The
new finite-stratum nuisance certificate protects unknown valid-subset means
and squared drifts, rather than an unprotected all-source average. At K2048,
weak bias and target heterogeneity .3, binary length/ordinary-target changes
1.736 to1.079 (hybrid) or1.066 (dispersion); no target floor gives .301.
Fitted K128 lengths improve1.110 to .872, but an equally improved full-target
certificate is shorter at .831. Useful fitted borrowing is not yet achieved.
All observed source-procedure coverages are1, so conservatism remains.
Artifacts: implementation/next_paper/v16/unified_inference_20261002_v1/.
Published to the dedicated RoCE-K Overleaf project at `1e15b0ec564b8ad94aa0c762f30faea53d578838`.
Main entry: `main.tex`; paired experiment report: `experiments_20261002b.tex`.
The current RoCE/ENAR project was not modified.

Use report/REPORT_zh.txt and compiled/concentration_note.pdf. The separate
paired_target_audit corrects spurious strict-zero floating-point diagnostic
flags; no intervals were changed. All previous artifacts remain intact.

NEXT-PAPER REPAIRS (2026-10-02): completed18 Gaussian settings x2000,
18 binary-score settings x1000,12 fresh saturated-nuisance settings x300,
plus a paired shared-target budget refinement. Reported total57600 repetitions;
the paired refinement is not counted again. All sites total1000 patients.
Additionally audited90 archived actual RoCE fits and40000 target-only draws.

The shorter protected fallback has a pointwise source-moment length cap.
At K2048,no shared target floor,all-valid midpoint RMSE changes .0074 to .0010;
at strong shifts .0076 to .0021. Fixed target error allocation remains preferred:
decreasing it can widen intervals with a target floor (.526 to .695 length ratio
in the weak-shift example). Binary Bernstein calibration reduces the illustrated
known-variance length ratio from2.056 to .924; empirical target variance gives
1.089, and source variance certification still costs1.737.

Fresh saturated fits expose the remaining barrier. At K128,oracle source drift
allowances average about .014, feasible budgets about .34, and paying shared
target uncertainty once improves this to about .28. The resulting fitted
interval is still about1.10--1.11 times an equally finite-sample-protected
full-target reference. This is not a completed efficient high-dimensional
procedure. Production RoCE and prior manuscript/results remain unchanged.
Artifacts: implementation/next_paper/v15/calibration_fallback_20261002_v1/.
Read report/REPORT_zh.txt or compiled/repair_note.pdf. Frozen sources, tests,
paired records, covariance conventions and the baseline audit are preserved.
The standalone report is published in the new RoCE-K Overleaf project at
c2c61db0b7d453616e2a0231ea4095acf347d6a3; entry experiments_20261002.tex.
Only six new report/assets files were added; existing manuscripts were retained.

NEXT-PAPER RESIDUAL EXPERIMENTS (2026-10-01): completed 114 reported settings,
138000 setting-specific repetitions, K=8--2048 and n=1000/site. This is an
oracle residual-score study, including Gaussian summaries and exact binary
patient-score simulations, not full high-dimensional nuisance refitting.
The target uses a joint three-component region; source sets use analytic
Fourier inversion, dispersion or coordinate counting. Independent confirmation
combines Fourier and dispersion with a prespecified split of error budgets.

At K=2048 without shared target noise, confirmation weak-boundary naive pooling
covers43.15%; combined certificates cover99.2%, with mean length16.23% of the
ordinary target interval. At strong shifts combined coverage is98.95%, length
14.30%. Dispersion is shorter at weak shifts; Fourier is shorter at strong
shifts. Known source identities remain an explicitly privileged oracle.
The protected procedures overcover; these results do not prove optimality.

Uncertainty costs remain material: the budget-corrected variance certificate
nearly removes the length gain in the illustrated shared-noise setting;
fully bounded binary-score protection can exceed target-only length.
Rare target fallbacks dominate midpoint MSE in several confirmation cells
(0.8% of all-valid hybrid repetitions contribute97.7% of squared error).
Fixed-probability target-width fallbacks also need control for any expected
length theorem. Actual25%/50%-valid boundary settings preserve multiple
source-compatible centers and show weaker borrowing gains.

All final reported sites respect1000 total patients. The variance panel was
rerun with250 source pilot+750 evaluation observations; its original six
extra-pilot cells are preserved but excluded from final summaries. Nine
implementation tests pass and984 method-setting rows are recomputed from
repeat records. Current-paper estimators and prior results are untouched.
Code: diagnosis/next_paper/residual_compatibility/.
Artifacts: implementation/next_paper/v14/residual_compatibility_20261001_v1/.
Read report_final/REPORT_zh.txt or compiled/experiment_note.pdf.

NEXT-PAPER METHOD DESIGN V2 (2026-10-01): published and remotely verified at
1e751ab6a2c5c930e53da194cb07fa88be861abe in the separate RoCE-K Overleaf
project 6abed6440abbc7c132b9611d. The primary direction is a fixed-dimensional
joint compatibility region from target means, source centers and unbiased
source dispersion, projected onto TATE. It is a conservative moment
relaxation of count compatibility, not Guo's full searching/sampling method.
One-dimensional borrowing is optional; the version-1 joint GLS manuscript
remains available as legacy_gls.tex. Current RoCE/ENAR defaults are unchanged.

There are 990 passing mathematical assertions, including 4096 exact discrete
score configurations and 100 joint-ellipsoid witnesses. The ideal Gaussian
radius diagnostic confirms the practical cost of Cantelli calibration; it is
not fitted coverage or expected length. Useful one-sided calibration, feasible
simultaneous nuisance budgets, estimated covariance and full growing-K
inference remain open. No new fitted simulation campaign or Slurm jobs were
started. Both manuscripts compile without warnings; new pages were inspected.
Canonical sources: diagnosis/next_paper/roce_k/.
Artifacts: implementation/next_paper/v13/method_design_overleaf_20261001_v2/.
The prior state and README are archived there; all old results are preserved.

LATEST ENAR RHC PUBLICATION (2026-09-28): published andverified at
c2d12b5cc2b204f16264e2e8033ea429f12d5bbd,baseadfac942 (latest user edits retained).
RHC now reports the audited outer-training preprocessing version:3.0158pp,
SE2.1643pp,CI[-1.2261,7.2577]pp,12.0614% lower reported SE than calibrated
Target-only. Body,weights,figure anddata agree. Outer preprocessing andshared
inner/CV transformation are explicitly described;SS/IVW andfull-sample DR
protocols are distinguished. Only six ENAR files changed. Body stays12 pages;
compilation has no warnings andfull RHC pages were visually checked.
Simulation,theory,main.tex andold assets remain intact.
Artifacts: implementation/r12/enar_rhc_outer_preprocessing_20260928_v1/README.md
Selected input: real_data/rhc/report_selection_private_min_outer_preprocessing_v1/

RHC OUTER-PREPROCESSING PAIR COMPLETED (2026-09-28): all10 tasks finished.
Private/min,all5039 patients,61 features,published fold memberships,Two-layer,
one-round,joint TATE,cutoff2,source radius3,target radius5 stay fixed.
Whole-cohort controls reproduce every displayed historical metric exactly;
full target/source/initial nuisance signatures also agree within1e-10.
New outer-training preprocessing: RoCE3.0158pp,SE2.1643pp,CI[-1.2261,7.2577]pp;
calibrated Target-only2.4002pp,SE2.4611pp. Reported SE reduction12.0614%.
Old RoCE2.9380pp,SE2.1660pp; point shift+0.07775pp. SS/IVW change slightly;
full-sample Federated-DR/Pooled-DR are unchanged. Actual cohort/fold identities
and all10 outer-training row sets are verified. No clinical bias/coverage claim.

Use run_rhc_tate_experiment(preprocessing="outer_fold") for this sensitivity;
"cohort" remains the historical default. Initial/calibration CV still reuse the
outer-training transformation; this is NOT fully nested preprocessing. Source,
target,OR andweight share one transform within each outer problem. SS/IVW retain
their own original folds andexclude the relevant site's validation rows.
Runtime collation is fixed to historical C.UTF-8; explicit category dictionaries
preserve the old reference levels. testthat otherwise forces C collation,which
explained the initial fingerprint failures.40 tests/1043 assertions passed.
Numerical R/C++ library is frozen at implementation/r12/rhc_outer_preprocessing_v3.
Workspace additionally declares withr as a test-only suggested dependency.

All10 jobs were initially submitted; still-pending jobs were cancelled before
serial execution on the already allocated4-CPU node. Per-task locks andexecution
receipts prevent duplicates. saffo-2tb had no eligible nodes because its two
nodes acl45/acl47 are excluded for documented I/O failures; routing was corrected.
This paired result is now published in ENAR as recorded above.
Report: real_data/rhc/preprocessing_pair_v1/reviews/complete_v1/README_zh.md
Detailed fold/model audit: real_data/rhc/preprocessing_pair_v1/reviews/audit_v1/

LATEST ENAR WORDING/LABEL REVISION (2026-09-28): published and verified at
26e2ac5af3f75812c24c19c62cc3f49aeb737f79, based on cec8f9d3.
The abstract now says "For each treatment arm" and specifies source-assisted
estimates of the target potential-outcome mean. ENAR text, caption and both
simulation figures use Scenario 1/2/3 instead of C1/C2/C3. Internal simulation
configuration IDs remain unchanged; figure provenance records the label mapping.
The regenerated numerical CSV is byte-identical. RHC, theory, methods, DGP and
all other text are unchanged. Compiled with no LaTeX warnings; new labels and
abstract page visually checked. Current PDF, scripts and receipts:
`implementation/r12/enar_scenario_labels_20260928_v1/README.md`.

ENAR FIGURE AND RESULT UPDATE (2026-09-28): published and remotely verified at
cec8f9d3feb0aa80805b4317893568eef1c71c34, on base87c7368b9d7ef9baa7459fc6a8e7c4b33efade12.
The manuscript now displays only p=200 simulation RMSE/coverage: C1--C3,
K=2/4/6/8, rho=0/.5/1/1.5/2, all200 paired replicates per cell. A shared
publication style replaces the Beamer style in both simulation and RHC figures:
Okabe-Ito colours, method-specific shapes/line types, light grids, consistent
six-method order and actual-page typography. Columns are scenarios; rows are K.
RHC now uses the audited Private/min source radius3, target radius5 result.
Body, confidence interval, source-weight table and discussion match that fit.
All four captions are shorter. Only calibrated Target-only is displayed.
Method/theory, DGP, main.tex and old assets are unchanged. Old p100/p200 figures
and all data remain archived. No fitting or simulation rerun was needed.
The360 p200 metric rows match their saved values within4.9e-15 export rounding;
RHC matches its audited M3 export within1e-13. The18-page ENAR PDF compiled
without warnings, and all three figure pages were visually checked.
Artifacts and reproduction scripts:
`implementation/r12/enar_clear_figures_p200_20260928_v1/README.md`.
The previous state is archived in that directory as previous_current_state.md.

RHC RECOMMENDATION IS NOW AUDITED: Private/min,original partition,all5039
patients,61 shared features,Two-layer/one-round,joint TATE,aggregation cutoff2,
source fitting/evaluation radius3 andtarget radius5. RoCE is2.938pp,SE2.166pp,
95% CI[-1.307,7.183]pp,versus calibrated Target-only2.385pp,SE2.456pp.
The SE reduction is11.8%;the interval still includes zero. All four matched
baselines andthe target anchor are numerically unchanged. A six-method forest
with the Beamer palette/grid/font/order has been exported andvisually checked.
The complete result package is
`real_data/rhc/final_assessment_private_min_radius3_v2/README_zh.md`;
the selected-fit export is `real_data/rhc/presentation_private_min_radius3_v1/`.
This recommendation is a finite-sample choice,not proof of actual clinical
bias or all calibration assumptions. ENAR publication now follows the update above.

The additional clean OR-branch check is COMPLETE under
`real_data/rhc/model_validation_or_branch_v1/`:200 repeats of O_moderate_mix.
The target/outcome/propensity functions are unchanged;source feature laws are
90% target plus10% original source. True merged weights are.071–9.321,inside
radius3,andlog-weight projection residuals>2.3 verify working weight-model
misspecification. The12 initial/calibrated population moment systems have
gradient residuals<1.36e-8,positive Hessians(minimum eigenvalue>.00417),and
clipping margins>.997. This fills the clean OR-correct calibration-prerequisite
check while retaining the original severe-shift boundary case.
All200 repeats andreviews2271639/2271738 completed. Radius3 bias=-.07365pp,
RMSE1.3111pp versus Target-only2.1326pp (38.5% reduction),coverage95.5%.
Target-only andoracle results match the original OR-law target sample in every
repeat. Its plot was rendered in a fresh shell andvisually reviewed. The final
version2 assessment includes this supplement;the original four-case version1
report remains intact. All1000 scenario/repeat pairs across five laws are
accounted for,andthe full result package passed its checksum audit.
A report-builder final-path variable shadowing error was caught before any
final directory was published,corrected,andthe successful build verified.
The failed temporary build is retained in `report_build_attempts/radius3_v2/`.

RHC DIAGNOSIS AND TRUNCATION: user authorized tighter truncation while retaining
Private,min,Two-layer,one-round,joint TATE,cutoff2. Source fitting/evaluation
radii2/3/4/5 were prespecified; target fitting/evaluation radius stays5.
16/16 discovery fits (fold0/101/202/303) andreview2263383 completed. Radius5
reproduces the saved reference; target and initial OR messages are fixed.
60/60 fresh fits (radii2/3/4,seeds401–420) andreview2263524 also completed,
with20 matching saved radius5 references reused. No patients were dropped.

Fresh source-radius results (TATE/SE ranges in percentage points):
r2:4.29–5.50 /1.68–1.84; r3:2.69–4.76 /2.11–2.37;
r4:3.25–8.21 /2.50–3.46; original r5:4.46–12.14 /2.84–6.26.
Largest individual variance shares:r2 6.4%,r3 16.2%,r4 43.4%,r5 76.0%.
Median mean-absolute SMD:r2 .173,r3 .129,r4 .127,r5 .144.
R3 is the current compromise candidate:19/20 partitions have lower SE than
calibrated target-only (median ratio.896),with improved ordinary balance.
R2 has smaller SE but worse balance; it must not be selected simply for
significance. R3 caps upper weights for4.5–5.0% andlower weights for18.5–24.2%
of source observations. Symmetric log clipping can introduce bias; no true
bias/coverage claim follows. The old radius5 presentation is archived; ENAR now reports the recommended radius3 package.

Original-min decomposition is checked:target2.3846pp + source increments
1.8209,1.3043,2.3850 =7.8954pp. Medicaid record423 hasY1,prediction.1452,
merged weight148.413 andeta.17527; its single held-out residual contributes
3.4366pp. Score-only neutralization with all fits/weights fixed gives4.4587pp,
not a deletion/refit estimate. Fixed-eta SE4.257 versus reported4.511 shows
extra eta variance is not the main issue. Training variance ratios.55–.63
fail to predict realized held-out concentration; optimizer KKT residuals<1e-10.
Evidence: `real_data/rhc/reviews/min_explanation_v1/`.
Discovery report: `real_data/rhc/source_truncation_v1/reviews/report_v1/`.
Confirmation report: `real_data/rhc/source_truncation_confirmation_v1/reviews/report_v1/`.
Code reuses the checked v2 estimator library. New source_truncation profile
and source-radii argument drive the existing worker; five launcher tests and
an actual target-isolation fit check pass. Upper/lower cap counts are separately
reported. The synthetic isolation fixture emitted glmnet small-class warnings;
all real-data fitting andpaired-reference review checks passed.

RHC CASE REFITS ARE COMPLETE: `real_data/rhc/case_refit_v1/` has all18 fits
and reviewed output `reviews/complete_v1/`. All nine full-cohort controls
reproduce saved estimates/SE, target/source coefficients and eta within1e-10.
The nine diagnostic removals preserve every remaining outer-fold membership
and original preprocessing, but fully refit/retune nuisances and eta. No primary
patient is removed. For published fold0, removing Medicaid423 changes TATE by
-4.57pp at radius5, versus-.55pp at radius3 and-.62pp at radius2. At fold412,
the corresponding changes are-6.23,-1.16,-.37pp. Target303 removal in fold419
changes TATE by-.47,-.82,-.53pp at radii5,3,2; target-only changes-1.12pp.
Thus smaller caps reduce the identified source-case instability but do not
eliminate all sensitivity. All recorded fitted-model convergence checks pass.

A frozen known-truth auxiliary study is prepared as `real_data/rhc/model_validation_v1/`:
four prescribed RHC-feature population laws,20 repeats per law,paired radii2/3/5
within each job. Population probability identities, fixed exact truth and oracle
expectation/variance checks pass. First four jobs2265385–2265388 finished:
O_case_mix,W_tail,W_departure complete;W_overlap fails an initial-weight CV.
The original76 remaining tasks were not released. Completed-fit resume of the
first task reproduces every saved method row without modifying saved fits.
The [design](repair/rhc_model_validation_design.md) explicitly distinguishes
these model-based diagnostics from actual RHC causal bias and formal coverage
validation. In the empirical-case-mix scenario, discrete site-specific feature
support creates severe true-weight clipping at radius2/3; it is deliberately
retained as a support stress test. Other scenarios have exact exponential
joint-weight models,with or without an invalid treated source. No manuscript
or clinical primary result is replaced while validation remains incomplete.

The W_overlap failure is reproduced with PN,CD andwith/without the CV certificate.
Its full-sample lambda_max is.1134,while CV training-fold bounds reach.1657.
An adequate upper endpoint restores convergence. `R/model_fitting.R` now retries
only an entirely failed initial-weight CV path,appending larger penalties based
on feature extrema while preserving every original candidate andthe CV folds.
Successful original paths/fixed lambda fits stay unchanged; `cv_grid_retry`
records the extension. No convergence standard or compiled solver changed.
Frozen implementation: `implementation/r12/initial_dr_cv_retry_v1/`;18 regression
tests/292 assertions pass,andthe exact failed RHC-model subset is recovered.
Full real-data parity job2265646 completed:both radii3/5 reproduce saved points,
SEs,model coefficients andeta within1e-10.

`real_data/rhc/model_validation_v2/` retains the exact population template and
replication seeds,with only this checked recovery implementation. First four
jobs2265695–2265698 all completed and passed the smoke review. Every previously
successful first-repeat fit (nine fits across three scenarios/radii) reproduces
the old point/SE with maximum difference0;models andeta also match within1e-10.
The formerly failing W_overlap completes all radii with recorded grid retries.
Release job2266296 passed the gates andsubmitted the remaining76 repeats:
2266759–2266835 excluding2266786. All80 scientific tasks are now submitted;
final review2266836 depends on the remaining jobs. Failed v1 results remain.
A consolidated Chinese
assessment andoriginal-partition table are in
`real_data/rhc/reviews/truncation_assessment_20260928_v1/`.

The full v2 batch finished with72 successes and8 failures,all in W_overlap
(repeat3/4/8/11/12/16/17/20). These are calibrated-weight CV grid failures.
For repeat12/source s3/treated/outer4,the original lambda_max=.03905091 is
below CV fold3's necessary constant-feature support floor=.04120191;extending
the grid restores convergence. `R/model_fitting.R` now uses one shared retry
helper for initial andcalibrated/refined density losses. See the
[CV recovery derivation](repair/density_ratio_cv_retry.md).
Frozen library `implementation/r12/density_cv_retry_v2/` passes20 tests/313
assertions andboth exact failure fixtures. Its v1 precursor is retained but
superseded only because its Rd help had duplicate details;v2 builds cleanly.

`real_data/rhc/model_validation_v3/` freshly refit all8 failures andthe four
successful first-repeat controls,without migrating caches. Jobs2267158–2267169
all completed. Fresh real-RHC radius3/5 parity job2267123 passed both radii.
Combined review2267170 passed all controls andclinical parity,then included
every original scenario/seed once,retaining both versions. All12 repeated
successful-control fits have point/SE differences exactly0.
The old full-batch review2266836 correctly refused to summarize missing repeats.

Interim complete-scenario results (20 repeats each) are saved under
`model_validation_v2/reviews/complete_scenarios_v1/`. W_overlap is wholly omitted
from that interim review,not averaged over its successful subset. At radius3,
RMSE in O_case_mix/W_tail/W_departure is1.767/1.401/1.487pp versus calibrated
Target-only2.173/1.367/1.367pp;coverage is18/20,19/20,19/20. These do not prove
uniform efficiency gains oradequate nominal coverage.

The W-scenario oracle's pilot SD is1.461pp despite exact population SD2.344pp.
An independent2000-draw oracle-only check reproduces all original oracle rows,
has SD2.391pp andcoverage94.5%;the O-case check has SD2.159pp versus exact2.137pp
andcoverage94.9%. This supports the generator/variance formulas andshows why
the small pilot's relative RMSE is imprecise;it does not establish RoCE gains.
Evidence: `model_validation_v2/reviews/oracle_sampling_check_v1/`.
The same frozen design was extended to200 per scenario using IDs21–200
(`--seed-start`),retaining the original20. All720 new repeats completed with
no unrecovered scientific failure. Combined review2270265 passed;all800
scenario/repeat pairs and4000 method rows reproduce their summaries under
`model_validation_mc200_extension_v1/combined_mc200/`. Oracle-only draws are
not merged into the fitted-method summaries.

Radius3 RMSE gains against calibrated Target-only are15.6%/41.1%/37.9%/31.2%
in O_case_mix/W_overlap/W_tail/W_departure;coverage is97%/97%/95.5%/96%.
In O_case_mix,bias is.171pp at radius3 versus.492pp at radius2. The paired MSE
differences against Target-only are negative in all four cells,with approximate
95% Monte Carlo intervals below0. Both the clinical andMC200 figures were
visually checked. No claim of actual RHC bias follows from these model-based
experiments.

A population coordinate audit also finds genuine tight-cap calibration
infeasibility in the severe case-mix law. At radius3,Medicaid treated income
reference-group h-weighted moment target=.031568 exceeds its maximum attainable
value=.015560. Radius2 violates more coordinates. Radius5 passes these necessary
coordinate checks,but that is not proof of joint feasibility orall assumptions.
Evidence: `model_validation_mc200_extension_v1/reviews/population_overlap_v1/`.
This boundary is retained in the recommendation,despite good finite-sample
coverage in the prescribed200-repeat experiment.

A research-only fold mass contraction was also checked in
`mass_stabilization_checks_v1/`: original radius5 point7.895pp/SE4.511pp becomes
4.022pp/first-order SE2.527pp. It does not beat the target anchor's2.456pp SE,
so it was not adopted orpromoted into the estimator;its inference remains
unvalidated beyond identity/anchor/influence-linearity checks. The original
calibration losses andjoint aggregation algorithm remain the recommended
procedure,with the explicit source radius3 setting.

ARCHIVED RHC PRESENTATION (radius5, superseded above): original min,Private target only,six displayed methods.
Retain only calibrated target-only and label it Target-only. It is2.38pp,
SE2.46pp,95% CI[-2.43,7.20]pp. Ordinary target-only remains an internal reference,
not a displayed comparator. Other methods are SS,IVW,Federated-DR,Pooled-DR,RoCE.
Source1se and rotated targets remain internal archives. No old result is deleted.

This selection is now published in ENAR at279146de0b490b2c2abfcaa021bce7ce2e635757
(base8f77d896),with verified remote receipt. RHC body,case-mix/weights table,TATE
figure and related abstract/introduction/discussion sentences were updated.
The figure matches the Beamer RHC style and six-method order; only the visual
style was reused,not old numbers. Latest source-index notation eta_{s_j} was
retained during the concurrent-edit merge. Methods,theory,simulations,main.tex,
Beamer.tex and old assets were preserved; no general concision occurred.

Private/min RoCE remains7.90pp,SE4.51pp,95% CI[-.95,16.74]pp on5039 patients,
61 covariates plus intercept andthree sources. The text acknowledges its CI
includes zero and attaches no source1se precision claim to min. Saved cohort,
arm contrast/variance,weights and all six-method CIs pass numeric checks. Final
30-page PDF compiles without warnings and figure/table pages were reviewed.
Artifacts: `implementation/r12/enar_rhc_private_min_20260927_v1/`.
Final PDF: rebased_beamer_v1/ENAR.pdf. Presentation package:
`real_data/rhc/presentation_private_min_calibrated_target_v1/`.
The earlier report_selection_private_min_v1 seven-method package is retained
as an internal audit source and no longer defines displayed target labels.

The joint TATE update to BOTH Overleaf main.tex and ENAR/ENAR.tex is published
and remotely verified at012b141729b4499849796d50ab520d94f3419b46, on latest
remote base c6a94124573cedba0354e9ea43f77a88c3dc2590. Concurrent abstract,
introduction and estimand edits were repeatedly rebased and retained. Main now
has arm-specific transport sets, two weight vectors jointly optimized for TATE,
full cross-arm covariance, signed unconstrained weights, the increment quadratic
and lambda_agg=1/cutoff. Its primary algorithm is Two-layer with explicit reused-
training-moment assumptions and a joint TATE oracle variance bound/proof. The
unchanged supplement is explicitly a supporting/optional nested and armwise-
constrained reference. ENAR adds the equivalent increment formulation and
clarifications while retaining its existing joint theorem.

Both PDFs compile; main has one pre-existing microtype footnote warning and no
new warnings; ENAR has no warnings. No calibration losses, simulation/RHC
results, preamble/layout or broader concision were changed. Main77pages;
final ENAR30pages reflects the user's concurrent introduction edit.12 saved
fits/120 objectives verify covariance, penalty, complete TATE and reported
variance; maximum objective difference1.82e-14. Artifacts:
`implementation/r12/joint_tate_manuscript_20260927_v1/` (validation and publish
receipts; final ENAR PDF in auto_sync_2/compiled_enar/ENAR.pdf).

ENAR simulation replacement is now published and remotely verified at
`ddb04632d9d5e2257983595b6eec96d8c2d268e7` (base246c5489). It uses the complete
p100/p200 MC200 main panels:24000 scenario-replicate datasets,840 metric rows,
seven methods and four updated vector plots. The existing DGP,method,theory,
RHC reference and manuscript layout are preserved. Results follow the existing
setting-to-RMSE/coverage narrative. Per the user's explicit latest direction,
no variance-conservatism discussion was added and the source-calibration
ablation paragraph was removed; no scientific result or old asset was deleted.
The user's new invalid textcolor annotation was minimally repaired while
retaining its pending-citation note. The31-page updated PDF has no LaTeX
warnings; all four integrated plot pages were visually checked. No broader
concision pass occurred. Artifacts: `r11/enar_simulation_revision_20260927_v1/`.

The September27 audit found no active FACE-HD jobs at its start. This was
not evidence that every experiment had completed: ten controllers had stopped
(four FAILED,six TIMEOUT). Rechecking all of their submitted workers found123
terminal I/O failures on acl45/acl47,including41 omitted by stale controller
ledgers in the earlier82-task audit. All123 are now resubmitted at their original
seeds after strict review and archival. All53367 checkpoints were read;53365
were readable,two empty files were preserved outside the resume inventory.
The recovery guard accepts only explicitly reviewed file/cache errors with
hash-pinned same-node/time corroboration; three focused tests passed.

Ten replacement controllers are RUNNING;first recovery verification has122
workers RUNNING and one valid committed result. Whole project snapshot:
132 RUNNING/142 PENDING. Scientific input hashes and worker caps/resources are
unchanged;66 routing overrides exclude acl45/acl47. Controllers stop after.5h
within1h allocations to reserve summary/requeue time. Recovery artifacts:
`r11/recovery_20260927_io_v1/`. These facts supersede earlier unresolved
simulation recovery/approval notes below. K8 supplementary method results are
still incomplete; its paired baseline campaign is now1200/1200 complete.

RHC target rotation is COMPLETE:48/48 new jobs and review2255741 succeeded.
Medicare1458, Private & Medicare1236, Medicaid647 each have3 rules x4 partitions
plus4 baselines; Private1698 results were reused. Source_all_1se in the discovery
Private & Medicare target has TATE5.19–6.64pp,SE1.87–1.97pp and gaps to the
calibrated target between-1.73 and-.64pp. This is a promising stability result,
not identified lower bias or a reason to choose the most significant target.

The all-target fresh-partition confirmation is COMPLETE:180/180 new fits and
review2258968 succeeded,plus60 reused Private fits. All240 method fits have
zero recorded nonconvergence or line-search failures. Full report:
`real_data/rhc/target_rotation_confirmation_v1/reviews/complete_v1/`.
For source_all_1se,SE is lower than calibrated target-only in20/20 partitions
at each target;median paired reductions are12.0% Private,41.3% Medicare,33.4%
Private & Medicare,43.8% Medicaid. These are reported-SE comparisons,not true
bias/MSE/coverage claims. Source-only1se largest-record variance shares reach
19.8%,26.5%,40.5%,20.2% respectively. The40.5% Private & Medicare outlier is
seed403,Private source record28,treated death,weight148.4; the earlier four-fold
maximum10.2% did not reveal it. Retain this finding in sensitivity reporting.

Earlier paper-readiness judgment (presentation scope superseded by latest user
direction above): usable as a real-data application,not as proof that bias is
small or all baselines are dominated. Keep the originally designated Private target as primary;
source1se is an exploratory tuning candidate. A focused high-influence case
sensitivity remains advisable before finalizing inference claims. No new jobs,
patient deletions,default changes or manuscript edits were made by this review.
Assessment: `target_rotation_confirmation_v1/reviews/paper_readiness_v1/`.

A report-only diagnostic defect was found for Medicaid's constant ortho=0:
using a1e-12 substitute target SD produced misleading SMD values. Corrected
SMD is NA for target-constant features, with raw weighted mean gaps separately
reported (worst ortho gap.00190 in the initial source_all_1se panel). Across
all48 matched method profiles every effect,SE,weight/ESS and individual influence
value is unchanged. The correction passed full saved-fit reconstruction at all
four targets. The active confirmation review is review_workflow_v2; old pending
review2257943 was canceled and replaced by2258968. No scientific job was canceled.
Corrected discovery report:
`real_data/rhc/target_rotation_v1/reviews/balance_corrected_v2/complete/`.

Rotating targets preserves each group's X/A/Y, feature basis/normalization,
and the original Private data hash. Interface112 assertions and Python4 tests
pass; generalized diagnostics reproduce all four historical Private tables
exactly. The existing driver now reads --target-site; diagnostics use actual
site mappings and minimum source-arm ESS fractions rather than assuming s3
is always Medicaid. Complete Private fits were reused, not selected anew.

Current bias-risk diagnosis: over20 confirmation partitions, source-only1se
minus calibrated target has gap1.15–3.64pp (median1.99pp); the Medicare-assisted
candidate gap is9.70–19.86pp (median13.33pp). Neither contrast is a bias estimate
because the target anchor is not known truth. Different targets imply different
estimands. Do not select the most significant target or infer smaller bias
from SE/ESS. Saved evidence: `target_rotation_v1/reviews/private_bias_risk_v1/`.

RHC source regularization diagnosis is COMPLETE: six rules times four common
outer partitions, 24/24 fits, on the historical 5039-patient / 61-feature cohort.
The new default reproduces the saved RHC primary with maximum numeric difference
zero. Target and initial OR messages are identical within each partition for
all source-only comparisons; selective final-stage isolation also passes.
Artifacts: `real_data/rhc/source_regularization_v1/reviews/report_v4/`.

Across these four partitions, original min has TATE4.03–7.90 percentage points,
SE2.83–4.51 and maximum individual variance share63.0%. The source_all_1se
candidate retains target min and has TATE4.67–6.28, SE2.18–2.36 and maximum
individual share7.8%. Medicaid treated ESS improves from25–82 to59–108.
Ordinary balance worsens: worst absolute SMD.708→.883 and median of mean
absolute SMD.144→.164. This is a precision/support tradeoff, not evidence of
lower causal bias or valid coverage. All24 fitted profiles have zero recorded
nonconvergence and final-fit line-search failures. Numerical CV certificate
exclusions are retained, not counted as valid candidate fits.

The prespecified fresh-partition confirmation is COMPLETE: seeds401–420,
min/source_all_1se/global_1se, 60/60 scientific jobs and review job2254248 all
succeeded. All input/workflow/library hashes match,80 fixed-stage checks have
maximum numeric difference0, and all60 profiles have zero recorded fit
nonconvergence or line-search failures. Artifacts:
`real_data/rhc/source_regularization_confirmation_v1/reviews/interpretation_v1/`.

For these20 new partitions, min TATE ranges4.46–12.14pp and SE2.84–6.26pp;
source_all_1se TATE4.33–6.55pp and SE2.02–2.46pp. Source-only1se reduces SE
and largest-record variance contribution in20/20 paired partitions. Median
SE ratio versus min is.630; versus calibrated target-only it reduces SE in
20/20,median reduction12%. Worst individual variance share improves76.0%→19.8%,
but mean absolute balance SMD worsens (median across partitions).144→.165.
The candidate is more stable, not proven less biased or best among all baselines.
Partition medians are descriptive, not a new pooled causal estimate. One cohort
cannot validate coverage/RMSE. No global default or ENAR RHC result was replaced.

The optional source CV controls share the existing calibration routine:
`calibration_control$source_lambda_rules` with initial_weight/weight/outcome.
Unspecified stages inherit nuisance_lambda_rule. v2 also forwards the same
controls through the existing source refit path. Installed-library regressions:
39 tests/1389 assertions; RHC interface59 assertions; manifest-grid3 tests,
all passed. v1 scientific fits and all old files remain immutable. Builds,
outputs and full checkpoints stay in scratch; only source/docs are in home.

Earlier covariate follow-up remains complete (36 named analyses): historical,
log_missing and grouped_log_missing primary results7.90/4.51,7.84/4.56 and
7.79/4.50 (estimate/SE in percentage points). These feature recodings did not
resolve concentration. Four baselines per representation and the repaired
ordinary-target result remain at their original paths. Completed report:
`real_data/rhc/reviews/covariates_complete_20260925_v1/`.
The replaced RHC status entries, including obsolete pending/approval notes,
are preserved in `source_regularization_v1/reviews/report_v4/previous_current_state.md`.

K8 fresh-seed study:74/1200 method and1135/1200 baseline results at23:52UTC.
Four controllers stopped after82 worker failures (60 K8,7 C4p200K2,12 half-
invalidp200K4,3 half-invalidp200K6); logs cluster on acl45/acl47 filesystem
I/O failures. Six additional older supplementary controllers timed out.
These have NOT been recovered in this evening audit. Evidence and exact
worker/controller identities are in `r11/current_paper_mc200_v1/reviews/
status_20260925_evening_v1/`. The checkpoint inventory contains22273 RDS files
(about69MB), including one zero-byte entry; read/validate before resuming,
quarantine only proven corrupt entries while preserving originals. Six workers
have explicit R file-connection error markers requiring the existing strict
hash-pinned reviewed_R_file_io path and corroborating same-node/time logs.
No failed seed may be discarded or replaced. Original p100/p200 main panels
remain complete and unchanged. The live project queue at the first audit was
7 RUNNING/128 PENDING; those counts include controllers and can change.

The user authorized further simulations and the current real-data application.
A separate K8 variance replication is now submitted under
`r11/k8_variance_validation_v1/`: p200,K8,C1-C3,rho0/1.5,fresh seeds201-400,
1200 method jobs plus1200 matching baseline jobs. Controllers2166888/2166890
are RUNNING. Original MC200 data and variance formulas are unchanged; review
new seeds separately before a clearly labeled pooled400 summary.

RHC now has12 submitted jobs under `real_data/rhc/current_two_layer_v1/`:
primary Two-layer/joint-TATE/cutoff2,4 baseline jobs,1se/radius12/source-standard/
ordinary-target sensitivities and three genuine outer-fold seeds. Primary
real-data radii remain5; no automatic fallback to a favorable profile. At the
first verification all12 jobs are PENDING (IDs in submission_receipts.json).
The RHC wrapper forwards current controls, reports mu1/mu0 weights separately,
and distinguishes ordinary/calibrated target references. A new fold_seed
explicitly changes partitions; the historical loader seed did not do so.
29 interface checks, existing RHC regressions and a fitted synthetic smoke
pass against the frozen library. R/real_data_rhc.R is a frozen interface overlay
for these jobs; the production estimation library is untouched. All workflow,
input and library hashes are pinned. The initial test helper's unintended
source build was moved to scratch; only v3 is the production interface gate.
No real-data result has yet been substituted into ENAR; concision remains
unauthorized. Preprocessing/timing/support caveats are recorded in the design.

The 2026-09-25 audit supersedes the historical progress counts below. Both
p100 and p200 main panels and their baselines are COMPLETE and paired:12000
prescribed data settings per dimension,200 seeds/cell. All frozen settings,
data hashes and ordinary-target references match. All source/target calibration
ablations and the covariate-shift2 study are also complete; the three new p200
calibration reports each validate600 pairs. C4 p100 has4000 complete pairs.
See `r11/current_paper_mc200_v1/reviews/status_20260925_v1/report.md`.

At p200, joint-TATE coverage under rho1/1.5/2 is93-98.5%, and RMSE is29.4-68.9%
lower than ordinary target-only. Weak rho.5 still undercovers, especially K2
(80.5-84.5%). K8 SE/empirical SD is1.16-1.21 at rho0 or strong departures.
Saved fixed-eta variance remains conservative (1.14-1.19), so the eta correction
is not the main cause. Source calibration alone still has small, uncertain
TATE MSE effects. Target calibration's p200 RMSE improvement is.6-1.3%, but
variance accuracy does not improve uniformly. Do not tune SEs to observed coverage.

Sixteen acl73 startup/I/O failures stopped ten supplementary controllers.
Every failed attempt was retained and resubmitted at its original seed; all64
campaign routes now exclude acl73. Scientific files and libraries are unchanged.
Ten replacement controllers were submitted after terminal-state and lock checks.
The first verification finds all26 replacement jobs PENDING, not completed;
project queue75 RUNNING/44 PENDING. Old NEEDS_REVIEW status files remain stale
until new controllers start. Recovery: `recovery_20260925_acl73_v1/`.
ENAR remains last verified at171c27e, with p100 figures; p200 reports are now
ready but have not been added to the manuscript. Concision remains unauthorized.
RHC has no current-method effect fit. Next: finish supplements, diagnose K8
variance and prepare the current-method real-data application.

The authorized ENAR update was published at171c27e on Overleaf main.
It retains Two-layer fitting, defines candidates by treatment arm, gives the
complete TATE and jointly optimizes two eta vectors using cross-arm covariance.
The theorem now states a joint TATE oracle variance bound against the calibrated
target-only contrast, with positive limiting curvature replacing the previously
unimplemented compact weight constraint. The text explains iid observations
versus correlated arm estimators and the correct-OR special case where the
population covariance is eta-invariant and separate/joint oracle optima agree.

The simulation section and two newly named figures now use the complete matched
p100 bounded_joint_v3 MC200 panel (12000 data settings,420 metric rows,7 methods).
Old figures and all old results remain; RHC is labeled historical. Twelve saved
fits/120 training objectives reproduce the full covariance, penalty and TATE
formulas within numerical precision. Four comparator routines match the frozen
library; their TATE construction is armwise fitting followed by a paired
contrast, including IVW's heterogeneity adjustment. The28-page final PDF has no
LaTeX warnings. Artifacts: `r11/enar_joint_tate_revision_20260924_v1/`.

The user explicitly withheld permission for concision until later approval.
No word/page reduction pass or manuscript font/spacing change was performed.
The earlier22-page draft and its concision budget below are historical; do not
execute that plan from these notes. See [covariance clarification](repair/joint_tate_covariance.md).

The user subsequently authorized the ENAR Method revision. It is now synced
to Overleaf main at615c8f8 (base d024679): only ENAR/ENAR.tex changed. The method
uses outer evaluation k1 and calibration k2; eta uses final calibrated fits on
the outer training sample. Necessary theory/discussion references were updated,
and training-moment control is an explicit condition, not inferred merely from
outer-fold independence. The user agreed that theory need not fix cutoff2;
the manuscript retains a growing cutoff and the existing simulations retain2.
The final22-page PDF compiled without warnings; the counted body/figure length
is19 pages against the verified2027 ENAR15-page requirement. The detailed
[concision plan](repair/enar_two_layer_revision.md) targets about4 pages of
substantive cuts and preserves the core proof chain. Its larger cuts and result
updates have not yet been executed. Artifacts and the verified publication
receipt are in `r11/enar_two_layer_revision_20260924_v1/`.

The historical morning audit on2026-09-24 counted10687/12000 p200 primary repeats:
K2/4/6/8 have3000/2996/2602/2089. K2 now has a complete3000-pair report under
`r11/current_paper_mc200_v1/reviews/main_p200_k2_paired_mc200_v1/`. Joint coverage
is94.5-97% under rho1/1.5/2 and80.5-84.5% at rho.5. At rho0, C1/C2/C3 coverage
is96.5/97/92%; C3's SE/SD is.953 and Wilson interval87.4-95.0%. Retain this
lower-coverage finding rather than claiming uniform nominal coverage.

A new overnight acl72 incident interrupted31 workers and stopped18 campaign
controllers. All31 failed jobs were on acl72, with startup or file I/O evidence.
The reviewed recovery excludes acl72 for all64 campaign routes and existing
pending project jobs. All31 workers were resubmitted at unchanged task IDs,
seeds and scientific settings; all18 replacement controllers were verified
RUNNING. Nine recovery workers already committed valid results at the first
check;22 were still running. Of1454 checked checkpoint RDS files,1450 were
readable; four zero-byte files were preserved separately and excluded from
resume. Every failed attempt and old receipt remains archived. Evidence and
current receipts: `r11/current_paper_mc200_v1/recovery_20260924_acl72_v1/`.
The narrowly scoped operator recovery guard passed two tests; frozen worker
code and estimation libraries are unchanged. Older September23 statements of
no controller failures are timestamped history, superseded by this incident.

The user also requested real-data follow-up. The RHC data-only preflight is
complete under scratch `real_data/rhc/preflight_20260923_v1/`: the historical
insurance cohort reproduces5039 rows,61 working features and one target plus
three sources with observed unequal sizes. Seven empirical source-arm support
flags and2651 missing urine-output entries require review. No new effect model
was fitted. The legacy RHC wrapper lacks explicit current two-layer/joint-TATE
controls and assumes one source-weight vector, so its interface and output
schema need alignment before a current-method application. See
[the review plan](repair/real_data_review.md). Simulations continue independently;
old real-data results are preserved and indexed by hash.

The new two-layer p100 main method is complete:12000 repeats,200/cell at
K2/4/6/8, C1–C3 and rho0/.5/1/1.5/2. The complete primary/target-reference
report is `r11/current_paper_mc200_v1/reviews/main_p100_methods_mc200_v1/`.
External baseline pairing is now complete for all K2/4/6/8:12000 pairs passed
data hashes, ordinary-target references and frozen-setting checks. See
`main_p100_all_k_paired_mc200_v1/`, including360 paired MSE contrasts against
six reference methods. The earlier K2/4/6 report and incomplete-report attempt
remain preserved as historical records.

The historical snapshot `reviews/status_20260923_v2/` counted8436/12000
committed p200 primary repeats. Controller1725832 had timed out after submitting
2688 tasks, all of which subsequently committed. Replacement1868949 was submitted
only after terminal/live-queue and lock checks; it was PENDING at this snapshot.
The remaining312 tasks will use existing frozen inputs. No workers were cancelled
or duplicated. Configuration, manifest and provenance hashes are unchanged;
the controller now reserves15min for summary/requeue. Original receipts/logs are
in that campaign's `controller_history/1725832/`. The other current controller
receipts were live or COMPLETED; status-file timestamps alone are not authoritative.

For p100 joint TATE, rho1/1.5/2 coverage is94–98.5% across K2/4/6/8. At
rho.5, coverage is82–84%,87–88.5%,89–91%,91.5–93% respectively. The one-invalid
source design dilutes invalid information as K increases; this is not growing-K
or fixed-invalid-fraction robustness. Against ordinary target-only, strong-shift
RMSE is27.5–67.5% lower in the fully paired K2/4/6/8 cells. At rho0, some pooling
baselines have slightly lower RMSE than the adaptive method; no universal
dominance or exact95% coverage is claimed. All unfavorable rows are retained.

The p100 Two-layer calibration review is now complete at200 seeds/cell. Main
source/target ablations and the covariate-shift2 source ablation each passed600
pairs, including unchanged target/source reference checks. Source calibration
changes joint-TATE RMSE by-0.24%/-0.58%/-0.03% in C1/C2/C3 under the main DGP;
shift2 changes are-0.54%/-0.52%/+0.30%. All total TATE paired MSE contrasts are
within two MCSEs. Target calibration has little TATE RMSE effect but moves
reported/empirical variance ratios from.848/.856/.934 to.880/.893/.966.
The paired variance-discrepancy changes are about3.0-3.9 jackknife SEs;
this is not exact variance calibration or universal improvement.

Saved-score follow-up passed1200 source pairs across both designs:2400 fitted
programs reproduce both arm means, TATE and both variance formulas within1e-12;
actual outer folds match and eta records exclude outer evaluation observations.
Candidate/eta exchanges show small effects and some opposing contributions.
The overview, uncertainty, two-page figure and complete unfavorable results are
in `r11/current_paper_mc200_v1/reviews/calibration_p100_overview_mc200_v1/`.
Five focused Python checks and three existing review tests pass; both main
ablations reproduce all nine corresponding full-method TATE metric rows.
No production estimator, frozen library, DGP or submitted scientific setting
changed. At that historical time p200 campaigns continued; see `reviews/status_20260923_v4/`.

The 2026-09-23 evening snapshot counts9662/12000 committed p200 primary
repeats (K2/4/6/8:2925/2813/2218/1706). Main p200 baseline campaigns are
complete for K2/4/6; K8 remains incomplete. All current controller receipts
are live or COMPLETED. Throughput is presently constrained by queue priority:
129 project jobs were running and837 pending at the snapshot. The K2 and K4
controllers are live PENDING after checkpoint/requeue, not failed; do not
replace them or duplicate their workers. No new complete p200 panel or RHC
effect fit has been reported.

## Current-paper priorities

The current priority is to complete the current paper's prescribed fitted-method experiments,
with paired calibration and two-/three-layer variance comparisons, while
keeping next-paper research separate. No production estimator or DGP was
replaced by the next-paper confidence-set prototypes.

The user's new decision is to use two layers for the current paper and reserve
three layers for next-paper growing-K work. A frozen MC200 expansion is now
submitted under `implementation/r11/current_paper_mc200_v1/`; coordinator1725777
releases the existing checkpointed controllers. Main: p100/200, K2/4/6/8,
C1–C3, rho0/.5/1/1.5/2,1000/site,200 repeats/cell. C4, half-invalid sources
and source/target calibration ablations have lower concurrency and execution
gates. The plan has87700 new workers and reuses7500 matched baseline repeats;
its total in-flight cap is2732. Submission is not completion. The older
three-layer point estimates and coverage summaries below are historical
comparisons, not completed results for this new two-layer panel.

New theory-motivated DGPs are authorized if needed to examine method mechanisms
and relative performance. They must have separate design/version records and
retain every old result. Current submitted scientific inputs remain frozen;
see [the MC200 design](repair/current_paper_mc200_design.md).

The separate covariate-shift2 calibration study is submitted under
`r11/theory_stress_v1/calibration_mc200/`, coordinator1732041. It adds3300
new workers and reuses300 baselines for p100/200, K2, C1–C3, rho0 and200
repeats per method/cell. Its80 in-flight workers follow primary execution
checks; both new plans together cap in-flight workers at2812. This is a named
stress design within the bounded family, not a replacement for the main data.

Additional K8/p10/rho1.5 checks passed the complete ten-fold,100-lambda full,
standard-source and ordinary-target programs. The v1 local test mistakenly
passed a target-calibration-only control to the ordinary-target branch; the
production worker was already correct. The failed attempt is preserved and
the corrected v2 passed. Saved-score reconstruction and shared-eta diagnostics
passed37 assertions in `r11/mechanism_checks_v1/`.

Population-gradient checks over10 prespecified C1–C3/shift/strength settings
passed16/24-node quadrature and finite-difference checks. Standard-source score
derivatives reach.0259 in C2 and.00853 in C3; calibrated derivatives are below
6.1e-14, and both programs are population-consistent in these transportable
union-model settings. This verifies the intended first-order sensitivity
mechanism, not a finite-sample MSE or coverage advantage. Evidence:
`r11/theory_stress_v1/population_gradient_checks/`.

`crossfit_layers=2L/3L` is implemented and checked (1445 R assertions), with the
advisor's outer-evaluation/secondary-calibration definition. All p10/p100/p200
layer campaigns are complete:2700 results,1350 matched data pairs. See
[layer design](repair/layer_validation.md) and scratch
`implementation/r9/layer_review_lowdim_complete_v1/`. The new MC200 profile
explicitly selects two layers. Small K and small nuisance dimension p are distinct settings.
The more detailed `r9/layer_variance_review_v1/` adds paired arm/covariance
decomposition, jackknife Monte Carlo uncertainty and a two-page PDF. All216
component groups reconstruct TATE variance within2.8e-19;42 Python workflow
tests pass. The complete low-dimensional panel does not show a clear
three-layer advantage, but50 repeats/cell do not establish equivalence.
The complete p100 comparison is `r9/layer_review_p100_complete_v1/`: maximum
relative RMSE difference1.91%, with all paired MSE differences within two MCSEs.
Some paired variance-accuracy differences are about2.2 jackknife SEs, so the
report does not claim that every variance difference is undetectable.
The final combined report is `r9/layer_review_p200_complete_v1/`, with both
p200 PDF pages visually checked. At p200, TATE RMSE differs by at most1.25%;
all nine paired MSE contrasts are within1.64 MCSEs. Coverage agrees in seven
cells and differs by one repeat in two. Important exceptions: Two-layer mu0
RMSE is2.64% lower in C3/K4 (3.70 paired MCSEs), while Three-layer variance
is closer to empirical variance there. TATE variance-discrepancy contrasts
reach2.47/2.65 jackknife SEs in C1/K2 and C3/K4. Neither layer uniformly
improves variance accuracy. Fifty repeats/cell do not establish equivalence.

The one/two-round communication comparison is now complete, independently of
the layer comparison: `r7/communication_ablation_final_review_v1/` has900 pairs
across p100/200, C1–C3 and rho0/1/2, with50 seeds/cell. All data hashes, ordinary
and calibrated target references, and checked scientific settings match.
Joint-TATE RMSE differs by at most0.477%; all paired MSE differences are within
1.81 MCSEs. Coverage matches in17/18 cells. This supports retaining one-round
as the main protocol; it does not prove equivalence or settle layer selection.
The new report reproduces54 previously saved two-round metric rows.

All 606 newly prescribed r8 calibration-ablation results and their reused
ordinary-target references have committed. The complete matched report is
`r8/calibration_ablation_final_review_v1/`: all four variants have 300 results,
50/cell at p100/200 and C1–C3. Calibration effects are mixed at this sample
size; do not claim universal improvement.

The density-CV certificate passed all24 full-path native comparisons. Its
explicit R control `nuisance_cv_certificate` is off by default. Integration v2
passed690 focused assertions and4122 full regression assertions; existing
133 small-class warnings remain. C2,p10 two-round ten-fold results and7280
nuisance records match the previous library, and a separate two-worker PSOCK
check verifies argument propagation and serial/parallel equality. The p200
complete-pipeline comparison has also passed:7280 nuisance records have zero
coefficient difference and all aggregation results agree;8896 failed fits were
certified. A p200 baseline repeat matches all five frozen baseline estimates
and SEs exactly. New repeat campaigns can opt in through the checked launcher's
`--nuisance-cv-certificate` flag;43 workflow tests passed. Existing campaign
libraries remain frozen; see [CV validation](repair/density_cv_certificate.md).

The 2026-09-21 late-evening scheduler audit found four controllers terminated
by a transient squeue timeout and 17 jobs held after launch failure. Four
replacement controllers and release of the original 15 workers/two controllers
were verified. The controller-only repair passed 41 workflow tests, retries
failed scheduler reads without submitting from stale state, and writes status
atomically. Worker code, seeds and scientific configurations remain frozen.
Evidence is under `r9/alignment_review_20260921_2330/`. Future progress checks
must verify scheduler state as well as status-file timestamps.

The p100 K4/6 panel is now complete: 1200 method/baseline pairs across
C1–C3 and rho0/.5/1/2. Report and figures:
`r7/source_count_p100_final_review_v1/`. Weak-deviation coverage improves with
K in this one-invalid-source design, but its treated-arm weight remains
nonzero; the fixed-invalid-fraction panels remain essential. The matching
p200 panel is `r7/source_count_p200_final_review_v1/`: all1200 method/baseline
pairs passed data/reference checks. Joint coverage is94–96% at rho1/2. At
rho.5 it is78–90% for K4 and88–92% for K6; treated-arm weights on the invalid
source remain.120–.143 and.094–.102 respectively. At rho2, all its treated-arm
fold weights are zero at tolerance1e-10; its valid control arm remains used.
The original source-count and layer controllers1647126/27 and final layer
report1675447 completed with exit0, verified in Slurm accounting.
The completed300-pair C4 report is `r7/both_misspecified_final_review_v1/`; at rho0,
joint coverage is84/88% for p100/200 with bias about.011. This is an explicit
both-models-wrong boundary, not an assumption-valid primary setting.

The completed p100 strong-departure extension is now combined with its first
50 seeds: `r7/p100_deviation_mc200_review_v1/` has 200/cell at rho1/2, with
joint coverage94.5–97%. Weak rho.5 retains50/cell. The first p200 departure
panel is complete under `r7/outcome_shift_p200_final_review_v1/`; its weak
coverage is68–76%. The p200 extension has now also finished: the combined
`r7/p200_deviation_mc200_review_v1/` has200 repeats/cell at rho1/2, with coverage
93.5–97.5%, RMSE.0209–.0239 and SE/SD.989–1.116. C3,rho2 changes from90%
in the first50 seeds to93.5% in all200; its95% Monte Carlo Wilson interval is
89.2–96.2%. This reduces the earlier concern but does not establish exact95%
coverage. All900 extension method/baseline data hashes and ordinary-target
references agree; the combined1350 repeats have no overlapping seeds.
With TWO
weakly incompatible sources among K4, the completed150-pair panel has
coverage52–64% and bias.032–.035 despite SE/SD above1. This boundary is in
`r7/half_invalid_k4_weak_final_review_v1/` and is not explained by merely
having too few sources.
The matching K6/three-invalid weak panel is also complete: coverage32–34%,
bias about.038 and SE/SD above1. Its150-pair report is
`r7/half_invalid_k6_weak_final_review_v1/`. Both fixed-fraction weak panels
are now finished and retained as boundary evidence.

Five slow existing campaigns have higher operational concurrency limits after
checking the 5000-job user association limit. The caps are384/384 for the
p200 strong-deviation/source-count campaigns,192 for Three-layer p200, and96
each for Three-layer p100/Two-layer p200. Current controller identities and
receipts are under `r9/concurrency_expansion_20260922_v2/`. Its v1 preparation
stopped at a normal controller requeue transition and took no action. All five
replacement controllers started; worker jobs and scientific inputs were
preserved. Increased in-flight capacity does not mean all queued jobs start
immediately; use the saved live scheduler snapshots for actual throughput.

## Calibration review follow-up

The [review audit](repair/calibration_evidence_status.md) distinguishes the
current bounded DGP from the older simulation section in the saved Overleaf
snapshot. Current fitted OR/weight bases are shared, C2 target propensity is
correctly specified, and the saved main-text orthogonality lemma exists.
The direct calibration comparison is now complete, with mixed finite-sample
effects. Full manuscript alignment and a claim of systematic calibration
benefit remain unresolved.

The new source/target factorial ablation is implemented through an explicit
source nuisance control, with complete inner validation. It passed 998 R
assertions and 29 workflow tests. The default C2 regression fixture matched
the previous estimates and parameters exactly. Artifacts are under scratch
`implementation/r8/`; see the [design](repair/calibration_ablation_design.md).
No completed Monte Carlo calibration advantage is claimed from these tests.

## Storage and reproducibility

Code and design notes remain in this repository. Results, logs, checkpoints,
compiled libraries, downloaded papers and figure artifacts belong under
`/scratch.global/zhan9381/FACE-HD/`. Do not overwrite frozen campaign workflows,
scientific configurations, committed repeats or earlier reports. Retries retain
original seeds and an explicit failure history. Never discard repeats based on
coverage, bias or RMSE. Match data hashes and target references before comparing
methods. The targeted ENAR Two-layer revision was explicitly authorized on2026-09-24;
other manuscript rewrites and submission to a competition are outside that scope.

## Current-paper configuration

The default DGP is `bounded_joint_v3`: bounded features, four active slopes,
fixed population truth, and 1,000 observations per site. See the
[DGP specification](repair/bounded_dgp.md). The new main profile uses ten folds,
two layers, one-round communication, fold-summed source/target calibration,
100 nuisance lambdas, proximal Newton tolerance 1e-10, and radius 12.
Coordinate descent remains available. Cutoff 2 corresponds to aggregation
lambda 0.5. `joint_tate` is the main analysis; `separate_arms` and `common_tate`
remain reported comparisons. Common-weight performance is not a release gate.

Cutoff 2 activates a soft penalty; it does not force an incompatible source's
weight to zero. The saved manuscript's separated-bias, growing-cutoff theorem
does not establish uniform selection under fixed cutoff 2 and local departures.
Fixed rho=0.5 is a finite-sample boundary setting, not an asymptotic root-N
sequence. The user's corrected word was “fault”, not a request to make rho=1
its default. Failure of consistent screening does not by itself establish
population nonidentifiability.

## Completed experiments and checks

All paths below are relative to scratch `implementation/r7/`.

| Study | Completed repeats | Scope |
|---|---:|---|
| `validation_highdim_mc200_cutoff2_v1` | 1,200 | p=100/200, C1–C3, rho=0, 200/cell |
| `baseline_mc200_v1` | 3,000 | p=10/20/50/100/200, C1–C3, 200/cell |
| `lowdim_deviation_mc50_cutoff2_v1` | 1,350 | p=10/20/50, C1–C3, rho=.5/1/2 |
| `lowdim_cutoff2_control_mc50_v1` | 450 | Matched rho=0 controls |
| `lowdim_deviation_baselines_mc50_v1` | 1,350 | Matched deviation baselines |
| `extended_validation_v1` | 1,650 method + 750 baseline | p=100 deviations, covariate shifts and ablations |

The 1,200 high-dimensional rho=0 method/baseline pairs passed data-hash and
ordinary-target estimate/SE/truth checks within 1e-12; 60 metric rows reproduced
prior summaries. Joint coverage is 92.5–97%, with RMSE about .019–.0207.
Some baseline SE/empirical-SD ratios are about .88 and coverage 91–92.5%; retain
these for variance diagnostics rather than rescaling reported SEs.

The completed p=100 outcome-shift panel has 50 repeats per cell. Joint coverage
is 96% in all C1–C3 cells at rho=1 and rho=2. At rho=.5 it is 72–74%, with bias
.0272–.0301. Matched covariate half/double-shift controls have 94–98% coverage.
See `p100_extended_final_review_v1/`; its paired comparisons include uncertainty
and are based on identical generated data. Fifty repeats give imprecise coverage
estimates and do not settle all inference questions.

The completed low-dimensional deviation report is
`lowdim_deviation_final_review_v1/`. Joint coverage is 70–78% at rho=.5,
92–96% at rho=1 and 90–94% at rho=2. Independent optimization checks on 270 saved
problems agree with the production weights within 6.56e-8, with KKT residuals
at most 1e-10. Weak-deviation weight retention is therefore not explained by
an observed optimizer failure in those checked problems.

## Completed historical current-paper expansions

All32 campaigns in the three plans below now have every prescribed completion
marker:14300 method/baseline results in total. This is a completion inventory;
it does not replace each report's data-hash, baseline and statistical checks.
The layer comparison is a separate2700-result plan and is also complete.

- `highdim_validation_v2`: 4,750 method and 2,950 baseline repeats. It extends
  primary p=100/200, rho=1/2 cells to 200 prescribed seeds and adds p=200 shifts,
  ablations, C4 boundary settings and identical-population controls. See the
  [high-dimensional design](repair/highdim_validation_design.md).
- `source_count_validation_v1`: 3,000 method and 3,000 baseline repeats for
  fixed K=4/6. It includes p=100/200 with one deviated source, plus p=100 controls
  with half the sources deviated. Each cell has 50 repeats.
- `source_count_weak_fraction_v1`: 300 method and 300 baseline repeats, adding
  rho=.5 to the p=100 half-invalid controls. See the
  [source-count design](repair/source_count_validation_design.md).

The latter two studies vary fixed K; they do not establish growing-K theory.
K counts sources, in addition to the target. Increasing K changes both total
sample size and source-tilt positions. Results must be stratified by K and the
number of deviated sources. All new settings have matching baselines and passed
nonduplication checks. Twenty-five workflow tests passed. First-seed release
gates check execution and integrity, never favorable statistical performance.

## Execution and recovery

One repeat remains one Slurm job, with checkpoints and requeue. Worker routing
is `preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb`.
The partition list is not a strict priority ordering. New workflows exclude
acl03/17/18/19/22/41/66/72/74/83/96 and cn1114 after documented startup/I/O failures.
Operational routing can be updated without rewriting scientific inputs.

The acl18 recovery preserved 23 pre-R failures: 15 were resubmitted after
verified termination; eight were moved using the same job IDs. All 23 committed,
and the affected baseline campaign completed 900/900. Audit records are in
`recovery_20260921_acl18_v1/`. Earlier restart-count inheritance was fixed at
submission; a real requeue probe verified exact checkpoint restoration and
continued duplicate-attempt rejection. Do not bypass ownership guards or retry
a live estimator because a monitoring command timed out.

## Next-paper research

Research code and derivations are under [next_paper](next_paper/). Versioned
artifacts are under scratch `implementation/next_paper/`.

The v1–v4 entries below are historical milestones, not live scheduler status.
The latest research handoff is `next_paper/v9/WORKSTATE.md`. Later work includes
Gaussian dispersion intervals, a finite Bernoulli likelihood lower bound and
a bounded-score reference; none establishes validity of the full fitted RoCE
procedure. The honest-candidate interface has initial C1/C2/C3 implementation
checks; its inference remains unvalidated. Completing current-paper experiments
has priority over extending these reference studies.

The next fitted-data diagnostic is complete under
`next_paper/v10/honest_population_validation_v1/`:90 prespecified diagnostics
at p4, K2/8/16, C1–C3, rho0/.5 and five training seeds/cell. Each site still
has1000 observations, split500/500 for training/evaluation. Only eight research
workers ran concurrently. The extended population-moment audit exactly
reproduced old packets and passed independent shifted-outcome checks; the K8
pilot also passed artifact-preserving resume and bias-projection checks.
This is actual nuisance fitting, not Gaussian reference sampling. It does not
yet validate intervals. The complete report is
`next_paper/v10/honest_population_final_review_v1/`, with visually checked
`figures_v2/`. Oracle-valid TATE drift reached.527 conditional SD in oneK16
training sample; estimated/actual fixed-weight variance ratios were.928–1.072.
Neither candidate construction uniformly reduced drift. These conditional
diagnostics do not establish unconditional bias or coverage.
See the [fitted diagnostic design](next_paper/honest_population_validation.md)
and [mean-versus-dispersion bias decomposition](next_paper/honest_bias_decomposition.md).

The [independent covariance-pilot diagnostic](next_paper/honest_covariance_pilot.md)
is now checked under `next_paper/v11/`:500training/250pilot/250evaluation per
site, preserving1000 total. Ninety focused assertions passed; all90 original
fits were reused and180 full-held-out packets reproduced exactly. Pilot-based
oracle-valid GLS has reported/actual conditional variance ratios.865–1.057
(median.978), but relative covariance errors reach.584. This establishes the
sample-independence interface, not feasible covariance bounds or interval
coverage. The current-paper estimator and defaults are unchanged.

The research interval now accepts externally justified `valid_score_bias`
bounds. Its default/zero-budget behavior matches the previous frozen reference
exactly.835 assertions and independent small-K geometry checks passed; see
[bias allowances](next_paper/valid_score_bias_allowance.md). The
[honest block-score CLT route](next_paper/honest_block_score_clt.md) connects
the site-score representation to a Gaussian reference under explicit conditional
moment assumptions.180 fitted covariance packets passed the factorization and
whitening checks. This does not establish feasible bias envelopes, same-sample
covariance calibration or full RoCE interval validity.

- v1: a proved marginal confidence-set reference under a specified lower bound
  on valid sources; 145 known-Gaussian settings with 1,000 repeats each.
- v2: unequal private source variances and a common target component, calibrated
  through conditional Poisson-binomial probabilities; 167 settings with 1,000
  repeats each. In 165 assumption-valid cells, known-factor coverage is
  94.7–98.6%. The two valid-count violations yield 36–45.5%. These are known-score,
  known-variance references, not fitted high-dimensional RoCE validation.
- v3: a Gaussian variance-region prototype, with 40 focused assertions and
  80 end-to-end random-input checks passed. In 79 draws the simultaneous region
  contained the true variances; all 79 satisfied the required domination check.
  Initial numerical integration failures and their exact inputs are retained.
  These are small implementation checks, not a coverage study or validation of
  fitted-score variance bounds. No current-paper estimator uses this prototype.
- v4: 232 known-variance procedure settings on 100 paired data settings, with
  1,000 repeats each, completed. Guaranteed valid fractions are 25%, 50% and
  75%; the study contrasts partial and all-valid quorum rules. All 1,624 metric
  rows reproduce saved intervals, and 132 cells exactly match paired baseline
  intervals. Known-correlation coverage with nonzero deviations ranges from
  95.1% to 98.4%; oracle-relative lengths can still increase with K. Exportable
  figures and full tables are in `next_paper/v4/`. A separate 4-by-1,000
  estimated-variance study was launched; its first n=100, K=4 cell was conservative
  (99.8% coverage) and does not establish efficiency.

The variance-region coverage argument and a conditional growing-K score
reduction are documented, but full uniform nuisance rates, valid variance
regions from fitted scores, arm-specific validity and efficiency remain open.
A sufficient score-level growth bound is not a completed RoCE theorem. Guo's
2018/2023 JRSSB and the 2025 JASA RIFL references were verified from primary
sources; RIFL already addresses selection-robust federated inference, so no
novelty claim is made without further comparison. The user handles the current
manuscript. Continue experiments and next-paper research without declaring the
broader research objective complete.

New research notes distinguish uniform coverage from local oracle adaptation:
[half-valid bound](next_paper/half_valid_lower_bound.md),
[majority adaptation bound](next_paper/majority_oracle_adaptation_bound.md),
[uniform score remainder](next_paper/uniform_calibrated_score_remainder.md), and
[published-paper scope audit](next_paper/literature_scope_audit.md).
The Gaussian lower bounds have independent numerical checks but still require
an observed-data embedding before being claimed for the full RoCE model.

# Weight-relearning bootstrap in simulations

## Current independent-sampling stage

Latest completed checkpoint: n=100. The science/audit/summary chain
`18353074` -> `18353177` -> `18353219` completed with all 101 final-batch
Slurm records COMPLETED 0:0. Independent review `18384884` verified all 100 raw
and 100 audit bundles, exact hashes/bindings, 8,400 method rows, 600 inference
rows and 7,200 numerical checks. Its statistical recalculation differs by at
most 4.885e-15. Maximum KKT/other errors are 1.976386e-7/4.996004e-16;
recorded hard nuisance and weight-bootstrap failures are zero. Saved-data tail
review 18385881 also completed 0:0: all 600 primary fits and 1,800 site cells
pass observation-level common-weight reconstruction, with maximum pseudo-value
error 3.553e-15 and variance error 2.168e-19. The largest tail remains
seed10029/rho1.5/target/row657 (centered -74.005842; 45.08% of total fitted
variance). The new-batch large positive case is seed10062/rho2.5/target/row74,
A=1/Y=1, centered +59.265627 (37.66% of total fitted variance). Its numerical
reconstruction passes; this alone does not diagnose its nuisance predictions.
There are 122 zero-SS site cells; no high-leverage observation was removed.

| Rho | RoCE RMSE | Empirical SD | Mean analytic SE | Mean relearn SE | Analytic coverage | Relearn coverage |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.024953 | 0.021979 | 0.021968 | 0.022856 | 93/100 | 95/100 |
| 0.5 | 0.025889 | 0.023656 | 0.021974 | 0.024648 | 87/100 | 91/100 |
| 1 | 0.032112 | 0.030474 | 0.023414 | 0.029497 | 82/100 | 88/100 |
| 1.5 | 0.033156 | 0.033320 | 0.025112 | 0.031413 | 83/100 | 92/100 |
| 2 | 0.031728 | 0.031189 | 0.025840 | 0.030740 | 90/100 | 93/100 |
| 2.5 | 0.030089 | 0.028559 | 0.026198 | 0.029740 | 93/100 | 94/100 |

Target-only RMSE is 0.034767, coverage 92/100, and mean error -0.011766
(MCSE 0.003288). All six RoCE RMSEs remain numerically lower, but paired MSE
improvements exceed two MCSEs only at rho=0, 0.5 and 2.5. At rho=1, RoCE's
mean error is +0.010575 (MCSE 0.003047), analytic SE/SD is 0.768, and relearn
SE/SD is 0.968. The diagnostic relearn coverage's exact pointwise 95% MC
interval is [0.799764, 0.936431]: improved variance scale does not validate
the interval. At rho=1.5 mean error is -0.000423 (MCSE 0.003332), yet analytic
coverage is 83/100 and SE/SD is 0.754. Bias and variance issues must be
distinguished rather than attributed to one SE multiplier. These multiple
pointwise/staged comparisons are not multiplicity-adjusted formal tests.
Neither interval policy is validated; do not promote the diagnostic CI or
launch final production under a claimed passed coverage gate. All adverse
seeds and original intervals remain unchanged.

Authoritative outputs are `summaries/n100/` and `review_n100_v1/` beneath the
pilot root. Summary checksum-manifest fingerprint:
`bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe`.
The independent `review_independent_n100.R` checks integrity and recalculates
statistics; it does not refit or add MC samples. The separate
`review_n100_influence_tails.R` and `tail_review_n100_v1/` provide the saved-data
tail reconstruction; neither is a proof of the asymptotic variance formula.

### Completed n=100 mechanism diagnostic

The subsequent saved-data mechanism diagnostic (18410137, 0:0) verifies all
600 exact target-plus-source decompositions. At rho=1, target bias -0.011766
plus mean s1 increment +0.022767 and s2 increment -0.000427 yields RoCE bias
+0.010575. At rho=1.5, bias is almost zero but mean analytic variance is only
57% of repeated-sampling variance. The full source, variance and diagnostic
studentization analysis is in `N100_COVERAGE_MECHANISM_DIAGNOSTIC.md`; its
truth-assisted thought experiments are not new intervals or validation gates.

### Historical n=50 gate and final-batch submission

The n=50 checkpoint had independent integrity, statistical
recalculation and primary-TATE influence-tail reviews completed. All tasks
1--50 / seeds 10001--10050 completed and passed numerical audits. The final
prespecified pilot prefix, tasks 51--100, is submitted as array `18353074`,
throttle 16, with 50 CPUs / 32 GiB per task. The CPU request changes only
scheduling: five CV threads, two source workers per arm and parallel arms stay
unchanged, while five positive-rho tasks can now run in one wave instead of
four plus one. No method, cutoff, seed or CI was changed. This is not validated
inference or overall experiment completion. See `checkpoint_n050_review.txt`
under the pilot root. The canonical production grid remains unsubmitted.

The current batch has an explicit postprocessing chain: audit array `18353177`
uses `aftercorr:18353074` for the same task IDs 51--100, so each audit can run
only after its corresponding simulation succeeds. Summary job `18353219`
uses `afterok:18353177` and runs the unchanged exact n=100 summarizer. Audits
use one CPU / 4 GiB each, throttle 4; the summary uses one CPU / 4 GiB. Both
require complete scientific inputs and refuse existing output directories.
Failed upstream jobs do not produce a passed checkpoint; inspect and diagnose
them rather than replacing seeds or duplicating the pending jobs.

The scheduling-only `run_independent_inference_pilot_audit_array.sh` passes
19 isolated wrapper tests, including strict array-ID/path mapping, incomplete
library/output/lock rejection and nonzero R exit propagation. A real audit of
the saved seed 10010 in an isolated temporary root reproduced all four
original audit payload hashes exactly, with no model refitting. Its first
sandbox attempt stopped at module logging permissions before R/output; the
permission-enabled retry passed. The frozen scientific package and nine-file
runner workflow are unchanged. See `checkpoint_n100_postprocessing_jobs.txt`
under the pilot root for IDs/fingerprints. Do not duplicate these automatic
audits/summaries while pending/running. After completion, independent integrity
and statistical review are still required; formal production is not auto-submitted.

### Complete n=50 descriptive checkpoint

The prior chain `18333300` -> `18333344` -> `18333380` completed 0:0. Exact
raw/audit/summary hashes and bindings pass independent review: 50 attempts,
4,200 method rows, 300 primary inference rows and 3,600 numerical checks.
Maximum KKT residual is 1.237465e-7; maximum other identity error is
4.996004e-16. Bootstrap and recorded hard nuisance failures remain zero.
Independent statistical recalculation matches the summary within 4.66e-15.

| Rho | RoCE RMSE | Empirical SD | Mean analytic SE | Mean relearn SE | Analytic coverage | Relearn coverage |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.025103 | 0.022288 | 0.021955 | 0.022686 | 46/50 | 47/50 |
| 0.5 | 0.025858 | 0.024500 | 0.022083 | 0.024987 | 43/50 | 46/50 |
| 1 | 0.030282 | 0.029829 | 0.023695 | 0.029963 | 42/50 | 46/50 |
| 1.5 | 0.031839 | 0.031649 | 0.025429 | 0.031283 | 44/50 | 46/50 |
| 2 | 0.031580 | 0.029811 | 0.026063 | 0.030337 | 45/50 | 45/50 |
| 2.5 | 0.030979 | 0.028110 | 0.026308 | 0.029135 | 45/50 | 46/50 |

Target-only RMSE is 0.033563 and coverage 47/50 at every rho. Its mean error
is -0.014942 with MCSE 0.004293, so the coverage question cannot be attributed
solely to weight learning. All six RoCE RMSE values are numerically lower;
paired MSE improvements exceed two MCSEs only at rho=0 and 0.5 (about 3.13
and 2.42). This is not a uniform or multiplicity-adjusted superiority claim.

Analytic coverage remains adverse at rho=0.5 and 1: exact pointwise 95% MC
intervals are [0.732604,0.941808] and [0.708874,0.928299], both with upper
endpoints below 0.95. Multiple rhos, intervals and staged prefixes are being
examined; these are not prespecified sequential/multiplicity-adjusted tests.
Relearned SE/empirical-SD ratios are now 0.988--1.036, but coverage is only
90--94%, and mean biases remain. Better variance scale does not establish a
valid interval. The original estimate/SE/CI and all adverse seeds are retained.

All 900 primary site cells also pass variance reconstruction (maximum error
2.17e-19). The new largest tail is seed10029/rho1.5/target/row657, centered
-74.005842, accounting for 67.56% of target-site SS and 45.08% of total fitted
variance. Its exact saved-data trace gives A=1,Y=0, propensity 0.02039244,
treated prediction 0.7056127 and control prediction 0.5240720. The target-only
TATE pseudo-outcome is -34.42014. Common fold-4 weights are s1=0 and
s2=0.288757692; target anchor 0.711242308 multiplies that large residual-driven
term, while the source target-side term is only +0.00213217 before scaling.
Multiplication by N/n_target=3 reproduces raw primary phi -73.436785765
exactly. This propensity exceeds the 0.01 floor; the case is not a clipping,
indexing or armwise-weight artifact. It is retained as a high-leverage case,
not used alone to explain population bias or to change clipping.

The largest source-site tails remain seed10016/s1 (+47.953252) and
seed10011/s2 (-38.812095). All 69 zero-SS cells are s1 cells with zero
contribution; their undefined concentration fractions are not invalid fits.
The only primary inference truncation remains the same seed10016 observation
across six paired rhos; safety clipping remains zero. Structured warning
capture stays incomplete for this v19 pilot.

Under the frozen protocol and completed implementation/statistical/tail
reviews, collecting the final n=100 pilot prefix is warranted, not promoting
the diagnostic CI. The launch adapter passed two positive mappings and five
rejected IDs; every new simulation/audit output/lock and n100 summary slot was
absent. The 40-to-50 CPU scheduling change preserves all scientific arguments;
current memory accounting has substantial headroom, but queue/total-runtime
improvement is not guaranteed and no dedicated p=100 four-versus-five-worker
replay was performed. Each positive task still computes its comparison methods;
only the RoCE TATE reuse path refits mu1/source1 and reuses mu0/source2.

### Completed tasks within the n=50 batch: seeds 10026, 10029 and 10036

These are individual execution/intermediate reviews, not an n=28 or n=36
statistical checkpoint. Scientific tasks `18333300_26`, `_29` and `_36`
completed 0:0 in 01:02:18, 00:54:34 and 00:53:20, respectively. Matching
automatic audits `18333344_26`, `_29` and `_36` completed 0:0. Raw and audit
payload hashes, raw-to-audit bindings and task/simulation/job mappings pass
review. Each seed has 84 method rows, six primary inference rows and 72
numerical checks. No implementation, hard nuisance or soft weight-bootstrap
failure is recorded; the soft diagnostic uses all 1,000 draws at each rho.

| Seed | Maximum normalized KKT residual | Maximum other numerical error | Comparison weights clipped / evaluated |
| ---: | ---: | ---: | ---: |
| 10026 | 3.650091e-8 | 4.996004e-16 | 10/2000 |
| 10029 | 1.237465e-7 | 4.163336e-16 | 10/2000 |
| 10036 | 3.508202e-8 | 4.163336e-16 | 51/2000 |

All three seeds' primary analytic and relearned diagnostic intervals contain
truth at all six rhos. Target-only and hard diagnostic intervals also contain
truth. These are single-seed observations, not validated coverage. Primary
inference-logit truncation and safety-clip counts are zero. Seed 10029's
rho=1.5 KKT residual is larger than the previous n=25 maximum but remains
below the frozen 1e-5 threshold; the other identities remain below 1e-10.
Review the KKT distribution at the complete n=50 checkpoint, without changing
the solver or excluding this seed. Its relearned/fixed SE ratios are slightly
below one at several rhos (minimum 0.957834); relearning weights need not
increase every conditional bootstrap variance.

Comparison clipping is distinct from primary TATE inference truncation.
The maximum source clipping fractions are 0.010, 0.008 and 0.049,
respectively. Logs for seeds 10029 and 10036 contain glmnet small-class
warnings; repeated console clipping messages are not additional independent
observations. Exact internal warning-call provenance remains unavailable,
so `warning_capture_complete=FALSE` and missing structured warning counts
must be retained. No bound, cutoff, CI or seed was changed during this review.

Tasks 27, 35 and 37 subsequently completed with their matching automatic
audits, also without implementation or bootstrap failures. Independent
raw/audit hash and binding review passes for each. The main review rechecked
all 72 numerical rows and six primary inference rows per seed: maximum KKT
residuals are 3.120885e-8, 3.934585e-8 and 3.936094e-8, with maximum other
identity errors 4.163336e-16, 3.885781e-16 and 4.440892e-16. All primary
analytic/relearned and target-only intervals contain truth; primary inference
truncation and safety-clip counters are zero. Comparison clipping is
33/2000, 17/2000 and 18/2000, respectively. Seeds 10027/10035 also have
small-class glmnet warnings in their logs. Seed 10037's low target-only point
(about 0.16215) is shifted upward by borrowing at every rho; this favorable
realization, like adverse earlier seeds, is retained without changing the
method. These per-seed reviews do not replace the pending complete n=50
statistical checkpoint.

### Adverse n=50-batch seed: 10030

Task `18333300_30` and its automatic audit `18333344_30` completed 0:0.
Raw/audit payload hashes and bindings pass independent review, all 72
numerical checks pass (maximum KKT 3.7230e-8), and no implementation,
soft weight-bootstrap or recorded hard nuisance failure is present. Primary
inference truncation and safety clipping are zero; comparison clipping is
also zero. Absence of stderr warnings does not retroactively establish
structured warning capture completeness.

Nevertheless, all six primary analytic intervals miss truth. Relearned
diagnostic intervals cover only rho=0.5, 1 and 1.5, and miss rho=0, 2 and 2.5.
The target-only estimate is 0.1606035 versus truth 0.2062656 (error -0.0456621),
although its wider interval contains truth. Exact fold-weight-preserving
source ablations were independently recomputed from the saved primary
Phase-3 components, including unequal target/source fold scaling:

| Rho | Source 1 increment | Source 2 increment | Final TATE error |
| ---: | ---: | ---: | ---: |
| 0 | +0.008385957 | -0.01987344 | -0.05714959 |
| 0.5 | +0.02485711 | -0.02213375 | -0.04293874 |
| 1 | +0.02237043 | -0.02488559 | -0.04817726 |
| 1.5 | +0.01440842 | -0.02676401 | -0.05801769 |
| 2 | +0.006943375 | -0.02766462 | -0.06638335 |
| 2.5 | 0 | -0.02683253 | -0.07249464 |

Zero weights reproduce target-only, and target plus these source increments
reproduces the final TATE within 2.78e-17. Source 2 aggravates the negative
target realization; source 1 usually offsets part of it, until its borrowing
vanishes at rho=2.5. These paired rho results reuse one dataset and are not
six independent observations demonstrating a population source bias.

The source-2 tail review identifies row 594, outer fold 4, A=1 and Y=0.
Its treated density logit is -2.595659 (within M=5), density ratio 13.40541,
outcome prediction 0.6639639 and residual -0.6639639, giving correction
-8.900711. The control-arm raw correction is zero for this treated row.
Under the actual common source-2 weight, its primary pseudo-value is
`3 * eta_s2 * (-8.900711)`, ranging from approximately -6.73 to -9.56.
It pulls the point downward and contributes to variance, but accounts for
only part of the source's net negative increment. No observation, cutoff,
source weight, nuisance parameter or CI was changed in response to this case.
It remains included for the complete n=50 statistical review.

Tasks 32, 38 and 39 also completed and passed the corresponding individual
hash/numerical reviews, with no implementation or bootstrap failures. Their
primary analytic/relearned and target-only intervals contain truth at every
rho. Comparison clipping is 31/2000, 14/2000 and 17/2000, respectively.
Seed 10032's log gives an unexpanded summary of 14 warnings; seeds 10038 and
10039 include small-class glmnet warnings. Missing internal warning context
is still retained explicitly, not interpreted as zero warnings.

### Further n=50-batch observations: seeds 10028 and 10040

Both seeds pass the raw/audit hash and binding reviews, all 72 numerical
checks, and the implementation/bootstrap/nuisance gates. Their adverse
interval outcomes are retained. Main-agent source ablations independently
reproduce target-only and the full estimator within 1.13e-17 for seed 10028
and 2.61e-17 for seed 10040, using actual Phase-3 fold weights and site scaling.

Seed 10028 has target-only 0.1328373, error -0.0734283; even target-only's
interval misses truth. Primary analytic intervals miss at every rho, while
relearned intervals cover only rho=0.5. Borrowing improves the point but does
not repair this large negative target realization: at rho=0, s1 contributes
+0.01939847 and s2 +0.00014817, leaving error -0.05388158. At rho>=1.5,
s1 borrowing is zero and s2 contributes +0.00068625, leaving -0.07274197.
The largest source tails are negative despite the small positive net source-2
increment; tail direction cannot substitute for the exact total increment.

Seed 10040 has target-only 0.1551044, error -0.0511612, but its interval
contains truth. The analytic interval covers only rho=0.5; relearned intervals
cover rho=0, 0.5 and 1, but not rho>=1.5. At rho=0, s1 contributes +0.01191396
and s2 -0.00662445, leaving error -0.04587167. At rho>=1.5, s1 borrowing is zero
and s2 contributes -0.01169325, worsening the error to -0.06285443. Neither
seed has an inference-logit truncation or safety clip. Comparison clipping is
29/2000 and 11/2000. Seed 10028's log contains an unexpanded 11-warning summary;
seed 10040 includes small-class warnings, still without exact fit attribution.

Seeds 10031, 10033, 10034 and 10041 also pass individual integrity/numerical
reviews, with all primary analytic/relearned and target-only intervals covering
truth. Seed 10031's relearned/fixed SE ratio reaches approximately 1.62 at
rho=1.5; this is retained as a conditional uncertainty diagnostic, not a
population variance estimate. No extra intermediate-prefix coverage summary
has been produced.

Later seeds 10042 and 10044 also pass their corresponding automatic and
independent reviews, with all primary intervals covering truth. Seed 10042 is
an important counterexample to expecting borrowing to help every realization:
target-only error is only -0.001580, but at rho=0, source increments -0.006745
(s1) and -0.021331 (s2) lead to final error -0.029656. At rho=2.5, source 2's
increment is -0.029815 and final error is -0.024235. This is a finite-sample
point-estimation loss, not proof that the nominally compatible source has a
population transport violation. Seed 10044's initially low target point
instead benefits from net borrowing. All these paired-rho realizations remain
included for the same complete n=50 checkpoint; none motivates deleting seeds
or tuning parameters after seeing their outcomes.

### Complete n=25 descriptive checkpoint

The prior computation array `18323368`, matching audit array `18323952` and
summary job `18323998` all completed 0:0. The exact prefix has 2,100 method rows,
150 primary inference rows and 1,800 numerical checks, with no missing,
replaced or failed seed. Raw/audit/summary payload hashes and reference bindings
pass independent review. The maximum normalized KKT residual is 4.65e-8;
recorded hard nuisance and soft bootstrap failures are zero. Independent
statistical recalculation agrees with the summary within 2.22e-15.

| Rho | RoCE RMSE | Empirical SD | Mean analytic SE | Mean relearn SE | Analytic coverage | Relearn coverage |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.022708 | 0.022396 | 0.021900 | 0.022596 | 24/25 | 24/25 |
| 0.5 | 0.029735 | 0.025804 | 0.021839 | 0.024535 | 20/25 | 21/25 |
| 1 | 0.032302 | 0.030332 | 0.023343 | 0.029640 | 21/25 | 23/25 |
| 1.5 | 0.032729 | 0.033403 | 0.024993 | 0.031070 | 23/25 | 24/25 |
| 2 | 0.031255 | 0.031323 | 0.025695 | 0.030150 | 24/25 | 24/25 |
| 2.5 | 0.029663 | 0.029026 | 0.025968 | 0.029021 | 24/25 | 24/25 |

Target-only RMSE is 0.034129 and coverage is 24/25 at every rho. All six RoCE
RMSE values are numerically lower, but only rho=0 has a paired MSE difference
exceeding two Monte Carlo SEs in magnitude (about 2.88); this is descriptive,
not a multiplicity-adjusted superiority test. Point-estimation performance and
interval calibration are different questions.

The rho=0.5 mean error is +0.015650 with MCSE 0.005161; at rho=1 it is
+0.012655 with MCSE 0.006066. The exact pointwise 95% MC interval for rho=0.5
analytic coverage 20/25 is [0.592963,0.931689], again with an upper endpoint
below 0.95. The relearned interval reaches only 21/25 there, with exact interval
[0.639172,0.954621]. Mean relearned SE/empirical-SD ratios are 0.930--1.009,
closer than analytic ratios 0.748--0.978, but this does not remove the positive
bias or all coverage deficits. Do not treat variance correction alone as a
solution. Multiple rhos, intervals and staged prefixes are being examined;
these are not preregistered sequential or multiplicity-adjusted tests.

The 450 primary-TATE site cells (25 seeds x six rhos x three sites) also pass
variance reconstruction, with maximum absolute error 1.08e-19. The largest
centered influence is seed10016/rho0/s1/row688 (A=1,Y=1): +47.9533, representing
60.6% of that site's centered sum of squares and 35.9% of total fitted variance.
Seed10011/rho1.5/s2/row498 (A=1,Y=0) is -38.8121, representing 39.24% of its
site sum of squares and 19.20% of total fitted variance. These finite,
high-leverage observations remain included; they do not alone establish bias
or an implementation error. Thirty-five source-1 site cells have zero
variance because borrowing is zero; their concentration fractions are
undefined, not nonfinite primary estimates.

During independent review, a proposed arm decomposition was corrected:
`arm_results$mu1/mu0$all_phi_agg` use their own armwise weights and cannot be
subtracted as a common-weight TATE decomposition. All final tail statistics
use the primary `fit$all_phi_agg`. Seed10012's known s2 row727 remains negative
(-16.9288/-20.2285 centered at rho=0.5/1), as independently read from that same
observation; a different row's armwise statistic must not replace it.

Under the unchanged protocol, the next exact prefix is n=50. Both execution
and statistical reviews support collecting that evidence, not adopting a new
CI or tuning cutoff/bounds. The stage-aware prelaunch test verified tasks 26
and 50 and rejected five invalid/out-of-batch IDs without fitting; all new
output/lock paths were absent. Only scheduling concurrency changes from 8 to
16. Original estimates, SEs, intervals, nuisance rules and all seeds are kept.

### Completed n=25 batch: seed 10012 review

Task 12 completed 0:0 in 00:52:29; its corresponding automatic audit
`18323952_12` completed 0:0 in 10 seconds. Raw/audit payload hashes and the
task/job/simulation mapping pass independent review. All 72 numerical checks
pass (maximum normalized KKT residual 2.78e-8), with no hard nuisance or
soft weight-bootstrap failure. This is a single completed task, not an n=12
checkpoint or a replacement for the complete n=25 prefix.

The rho=0.5 analytic and relearned diagnostic intervals both miss truth;
at rho=1 only the analytic interval misses. Target-only estimates 0.237829976
against truth 0.206265564, an error of +0.031564412, but its interval covers.
An exact foldwise reconstruction, respecting the actual target/source fold
weights and sample-size scaling, gives these source increments relative to
target-only:

| Rho | Source 1 increment | Source 2 increment | Net borrowing | Final point error |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | +0.032935837 | -0.003535084 | +0.029400753 | +0.060965165 |
| 1 | +0.025356461 | -0.004328536 | +0.021027925 | +0.052592337 |

Zero source weights reproduce target-only exactly; the full-minus-target
reconstruction error is at most 2.78e-17. Thus source 1 compounds the positive
target realization, whereas source 2 offsets part of it. Do not replace the
foldwise calculation with average reported weights times pooled estimates.

Source 2 row 727 (outer fold 1, A=0, Y=1) is a legitimate high-influence case:
density logit -3.71797 gives ratio 41.1807 within M=5, outcome prediction is
0.38313, residual 0.61687 and the signed TATE correction is -25.40302.
Its centered TATE influence is about -16.93/-20.23 at rho=0.5/1, accounting for
23.5%/24.0% of that source's variance. Its signed contributions to the final
point estimate are -0.005646679/-0.006747554. It enlarges the variance/tail
and pulls TATE downward; it does not explain the positive point error.
No source/arm/index mismatch or primary inference clipping was found, and no
bound, cutoff or observation was changed in response to this case.

Comparison DR diagnostics report 14/2000 weights clipped (0.7% overall,
maximum source fraction 1.4%); the log specifies 13 below and one above the
comparison bounds. It also says there were 16 warnings without storing their
full condition records. Preserve that limitation rather than inventing exact
call attribution or treating warning capture as complete.

### Completed n=25 batch: seed 10016 review

Task 16 completed 0:0 in 00:49:15 and its automatic audit `18323952_16`
completed 0:0 in 9 seconds. Independent raw/audit hashes and task/job mapping
pass; maximum normalized KKT residual is 2.76e-8, with no implementation,
hard nuisance or soft weight-bootstrap failure. Both primary intervals and
target-only cover truth at all six rhos in this seed, which is not a population
coverage validation. Relearned/fixed SE ratios reach 1.403/1.449 at rho=1/1.5.

Each rho explicitly records one inference logit truncation out of 2,000
(fraction 0.0005), maximum absolute raw logit 5.931321 and zero inference
safety clips. This is the prespecified M=5 truncation, not an unrecorded
numerical replacement or an automatic hard failure; do not count the paired
rho records as independent truncation events. Comparison DR separately
records 14/2000 lower-bound clips (0.7% overall, 1.4% maximum source fraction).
Small-class glmnet warnings remain in the logs without exact call provenance;
warning capture remains incomplete. No bound or nuisance rule was changed.

### Other completed tasks in the n=25 batch

Tasks 11, 13, 14, 15, 17 and 18 also completed 0:0 and passed their corresponding
automatic and independent per-seed audits. Each has 84 method rows, six
inference rows, valid raw/audit hashes and matching task/job provenance. The
recorded hard nuisance and soft bootstrap failure counts are zero. Both primary
intervals and target-only cover truth at every rho within these seeds; these
per-seed observations do not substitute for the planned n=25 statistics.

| Seed | Maximum audit KKT residual | Comparison DR clips / 2000 |
| ---: | ---: | ---: |
| 10011 | 3.27864e-8 | 18 |
| 10013 | 4.07279e-8 | 47 |
| 10014 | 3.54442e-8 | 39 |
| 10015 | 3.01627e-8 | 5 |
| 10017 | 3.33621e-8 | 10 |
| 10018 | 3.91011e-8 | 17 |

Seed 10011 has a source-2 maximum absolute centered influence of 38.81 at
rho=1.5, attained by a negative TATE value at row 498 (99th percentile of the
absolute centered influence about 5.715, site SD 1.960). Its source-2 variance contribution
is about 0.0004265 versus the target's 0.0004447. The actual standard error
0.029524 is finite and the variance reconstruction passes. The exact net
source-2 point increment is only +0.002465: a large tail/variance contribution
does not by itself establish a large net point shift. Target-only is below
truth in this seed, and realized source borrowing moves it toward truth.
Retain this tail for the complete checkpoint review; do not tune bounds to
remove it.

Seed 10013 logs glmnet error code -92: the 92nd lambda did not converge within
maxit=1,000,000, and solutions for larger lambdas were returned. This is a
path-truncation warning, not direct proof of a failed selected model, but the
saved objects do not retain a complete glmnet path/jerr/call binding. Therefore
the exact method/site/rho and selected-lambda distance from the truncated tail
cannot be established read-only. Do not interpret zero recorded hard counters
as proof that every glmnet CV path converged. A subsequent bounded replay is
described below; no fitting-rule change was made. The comparison clipping log specifies
46 lower and one upper clip, distinct from primary TATE inference.

Small-class glmnet warnings in seeds 10011, 10013, 10015 and 10018 also retain
incomplete call provenance. None of these cases was removed or relabeled as
warning-free merely because the saved estimates and audited identities pass.

### Seed 10013 target-only CV replay

Diagnostic job `18331540` completed 0:0 in 81 seconds. Its 25 PS and 50 OR
target-only complement-CV calls reproduced every saved prediction exactly
(maximum absolute error zero); the shared PS predictions also match the other
arm. All five output payload hashes pass. The replay produced zero warning
events, so it did not reproduce the original -92 warning within these calls.
This narrows the investigation but does not identify the original warning's
method/site/call among other fitting paths or establish that all original CV
paths converged. Original CV partitions were not saved, so the replay explicitly
does not claim a direct saved-partition identity check.

The standalone replay and its 46 no-fit test assertions enforce the original
NULL-group `nfolds` branch, seeds, family, penalty count, max iterations and
clipping rules. It preserves actual replay fold IDs, CV path/full-model jerr
metadata, warning events and warning-producing CV objects; no source model or
aggregation was rerun. Outputs are under
`independent_target_cv_replays_v19/seed_010013_v1` in the results base, separate
from the immutable simulation bundle. A warning lambda index is not assumed
to be an index in the master CV grid, since glmnet can generate separate fold
paths when its lambda argument is NULL.

Postprocessing job `18322944`, submitted with `afterok:18304069`, completed
0:0 in 17 seconds after task 10 completed 0:0 in 01:00:47. It ran the existing
seed-10010 audit and then the exact n=10 summarizer, with stop-on-error behavior,
one CPU and 4 GiB. No model was refitted. Both immutable outputs and their
payload hashes passed subsequent independent review; do not rerun into those
same output directories. The submission record remains
`checkpoint_n010_postprocessing_job.txt` under the pilot root.

### Complete n=10 descriptive checkpoint

All ten exact raw bundles and bound audits pass hashes/provenance: 840 method
rows, 60 primary inference rows and 720 numerical checks, with no missing,
replaced or failed seed. Maximum normalized KKT residual is 4.65e-8; hard
nuisance and soft weight-bootstrap failure counts are zero. Independent
recalculation agrees with all numerical summary fields within 4.58e-16.

| Rho | RoCE RMSE | Empirical SD | Mean analytic SE | Mean relearn SE | Analytic coverage | Relearn coverage |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.027070 | 0.027507 | 0.021499 | 0.022183 | 9/10 | 9/10 |
| 0.5 | 0.031743 | 0.029050 | 0.021592 | 0.024479 | 7/10 | 8/10 |
| 1 | 0.037692 | 0.036734 | 0.023063 | 0.029604 | 8/10 | 8/10 |
| 1.5 | 0.039664 | 0.041709 | 0.024634 | 0.031125 | 8/10 | 9/10 |
| 2 | 0.038152 | 0.040068 | 0.025382 | 0.030582 | 9/10 | 9/10 |
| 2.5 | 0.035827 | 0.037242 | 0.025660 | 0.029195 | 9/10 | 9/10 |

Target-only RMSE is 0.035300 and coverage is 9/10 at every rho. RoCE RMSE is
numerically higher in four settings; all paired MSE differences remain less
than twice their Monte Carlo SE. Analytic-SE/empirical-SD ratios are
0.591--0.782, and relearned-SE ratios 0.746--0.843. These are adverse signals
that must not be dismissed simply because n is small or because software
checks pass. The exact pointwise 95% Monte Carlo interval for the rho=0.5
analytic count of 7/10 is [0.347547,0.933260], whose upper endpoint is below
0.95. This is not a preregistered sequential/multiplicity-adjusted hypothesis
test: multiple rhos, intervals and checkpoint prefixes are being examined.
The 8/10 and 9/10 exact intervals are [0.443905,0.974789] and
[0.554984,0.997471]. The direct paired relearn-minus-analytic coverage changes
are +0.1 with MCSE 0.1 at rho=0.5 and 1.5, zero otherwise.

Under the fixed protocol, these statistical findings warrant more independent
seeds, not post-hoc tuning or favorable-seed selection. The next review is the
complete n=25 prefix, with the same method, nuisance rule, cutoff, bootstrap
counts and CI definitions. Only array concurrency increases from 4 to 8.
The stage-aware prelaunch test checked task mappings 11 and 25 and rejected
five invalid/out-of-batch IDs without fitting; all 15 new output/lock paths
were absent before submission. The scheduling adapter itself and all nine
scientific dependencies are unchanged. Inference remains unvalidated.

Seed 10010 retained comparison DR lower clipping of 27/2000 (1.35% overall,
2.7% maximum source fraction), distinct from primary TATE inference clipping.
Its log also contains small-class glmnet warnings. Every saved-data outer-
training arm/outcome minority count is at least 13; therefore those warnings
cannot be attributed to a full outer-training stratum with fewer than eight
observations. A smaller nested CV subset or another internal fitting call is
plausible, but exact caller/fold attribution is unavailable. Warning capture
remains explicitly incomplete, not zero. Both per-seed intervals contain
truth at every rho; that does not override the n=10 adverse findings.

### Per-seed audit history for tasks 6--9

In the n=10 batch, task 6 / seed 10006 completed 0:0 in 01:04:11 and passed
its 72-check real numerical audit and payload checksums. Maximum KKT residual
is 4.20e-8; implementation, hard nuisance and bootstrap failure counts are
zero. Its two intervals contain truth at every rho, which is a per-seed result,
not a new coverage validation. Comparison DR clipping is 11/2000 (0.55%
overall, 1.1% maximum within a source); small-class warnings are retained
without a fabricated exact call attribution. Task 10 started through the
existing array after task 6 released a slot; no task was resubmitted.

Task 7 / seed 10007 completed 0:0 in 01:09:26 and passed its 72-check audit
and payload hashes (maximum KKT residual 2.44e-8). It is nevertheless an
important adverse statistical case: both analytic and diagnostic relearned
intervals miss truth at every nonzero rho. Target-only already estimates
0.2670595 against truth 0.2062656 (error +0.060794), and its own CI misses.
The exact foldwise source increments relative to target-only are:

| Rho | Changed source s1 | Stable source s2 | Net increment |
| ---: | ---: | ---: | ---: |
| 0 | -0.020293 | -0.008898 | -0.029192 |
| 0.5 | +0.007199 | -0.006409 | +0.000790 |
| 1 | +0.024188 | -0.007470 | +0.016718 |
| 1.5 | +0.024456 | -0.008593 | +0.015864 |
| 2 | +0.020922 | -0.009225 | +0.011698 |
| 2.5 | +0.014226 | -0.009950 | +0.004277 |

Entries are rounded; the unrounded foldwise reconstruction error is at most
2.8e-17. At rho=0.5 the source contributions almost cancel, so noncoverage
primarily reflects the high target realization. At rho=1--2.5, residual s1
borrowing further raises the estimate while s2 offsets part of that increase.
The s1 mean fold weights at rho=0,0.5,1,1.5,2,2.5 are approximately
0.247,0.335,0.282,0.199,0.143,0.086. Wald statistics and penalties increase
with positive rho, but finite borrowing remains. Hard screening at rho>=1
equals target-only and also misses; it cannot be selected post hoc to repair
this realization. No hard nuisance, optimizer or bootstrap failure was found;
maximum centered primary influence magnitudes remain below 11.7. Comparison-
only DR lower clipping is 29/2000 (maximum source fraction 2.9%); saved output
does not identify the exact source/call. This seed is retained unchanged.

Tasks 8 and 9 / seeds 10008 and 10009 completed 0:0 in 01:14:17 and 01:10:50.
Both pass the 72-check audit with zero implementation failures; maximum KKT
residuals are 3.79e-8 and 3.13e-8. Seed 10008 misses with both intervals at
rho=0 (upper bounds 0.205685 and 0.205897, below truth 0.206266), while both
cover at positive rho. Seed 10009 misses with the analytic interval at rho=0.5
but covers with the relearned diagnostic interval; both cover at other rhos.
These paired, per-seed observations are not independent coverage estimates
across rho and do not validate the diagnostic SE. All original rows and CI
definitions are preserved.

Independent read-only review also verified both raw and audit payload hashes.
For seed 10008 at rho=0, target-only is 0.160367 with SE 0.03027: its wider
interval covers, while RoCE is only 0.00196 higher and its narrower intervals
miss slightly. Primary fold/site influence checks show no abnormal
concentration (largest single-point site-variance share 6.7%). For seed 10009
at rho=0.5, RoCE is 0.03460 above target-only (0.215279); the diagnostic SE
increase covers truth but remains a diagnostic, not a selected replacement.
Comparison DR lower clipping is respectively 37/2000 and 11/2000, with maximum
source fractions 3.7% and 1.1%; this is not primary TATE inference clipping.
All hard nuisance and soft weight-bootstrap failure counts are zero.
Warnings remain in scheduler logs without complete structured cross-worker
capture (`warning_capture_complete=FALSE`, `warning_count=NA`). Preserve that
limitation and address final-production warning provenance before claiming
complete warning auditing; do not mutate this frozen pilot to add it mid-run.

### Historical completion of the n=5 prefix

Within the earlier tasks-2--5 batch, task 3 / seed 10003 completed 0:0 in
00:57:56 and passed its 72-check real numerical audit and payload checksums.
The maximum normalized KKT residual is 4.09e-8, implementation failures are
zero and all 84 estimates/SEs are finite. Task 2 / seed 10002 then completed
0:0 in 01:21:35 and passed the same audit/checksums, with maximum KKT residual
3.53e-8 and no implementation failures. Completed, audited tasks are now 1--3;
task 4 / seed 10004 then completed 0:0 in 01:14:46 and passed its 72-check
numerical audit and payload hashes (maximum KKT residual 2.79e-8, no
implementation failures). Task 5 / seed 10005 subsequently completed 0:0 in
01:00:37 and passed its numerical audit/checksums, with maximum KKT residual
3.21e-8 and no implementation failures. The complete n=5 checkpoint was then
computed; no smaller subset was substituted for the predeclared prefix.

### Complete n=5 descriptive checkpoint

All five exact raw bundles and five bound seed audits pass hashes/provenance.
There are 420 result rows, 30 primary inference rows and 360 numerical checks,
with no missing or replaced seed. The maximum KKT residual is 4.65e-8; all
hard nuisance and weight-bootstrap failure counts are zero. The fixed TATE
truth remains 0.2062655638596 across seeds and rhos.

| Rho | RoCE RMSE | Empirical SD | Mean analytic SE | Mean relearn SE | Analytic coverage | Relearn coverage |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.026379 | 0.023635 | 0.021519 | 0.022302 | 5/5 | 5/5 |
| 0.5 | 0.023888 | 0.025439 | 0.021571 | 0.024504 | 4/5 | 4/5 |
| 1 | 0.030163 | 0.033131 | 0.023030 | 0.030325 | 4/5 | 4/5 |
| 1.5 | 0.034402 | 0.037384 | 0.024757 | 0.032166 | 4/5 | 5/5 |
| 2 | 0.033035 | 0.032164 | 0.025688 | 0.031376 | 5/5 | 5/5 |
| 2.5 | 0.031755 | 0.028521 | 0.025890 | 0.029263 | 5/5 | 5/5 |

Target-only RMSE is 0.033572 at every rho. Every paired MSE difference has
Monte Carlo SE at least as large as its absolute point estimate, so this small
sample does not establish superiority or inferiority. A 4/5 coverage count
has an exact 95% Monte Carlo interval about [0.284,0.995], and 5/5 about
[0.478,1]; neither validates nominal 95% coverage. The direct paired change
in coverage is +0.2 with MCSE 0.2 at rho=1.5 and zero in the other settings.
The seed-10004 noncoverage cases remain included without parameter changes.

The mean analytic-SE/empirical-SD ratios span 0.662--0.910 and the diagnostic
relearned-SE ratios 0.860--1.026. These are descriptive, not a validated
variance correction. The summary preserves `inference_validated=FALSE` and
`statistical_review_required=TRUE`.

Comparison DR lower-clip counts per seed are 23,17,11,13,45 out of 2,000 source
weights each; the largest single-source fraction is 4.5% in seed 10005.
Repeated method/rho rows are not counted as distinct weight fits. Small-class
warnings remain in logs without complete call binding; warning capture stays
explicitly incomplete rather than being recorded as zero.

The n=10 batch used the same method, nuisance rule, cutoff, bootstrap counts
and CI definitions. Only scheduling concurrency changed from 2 to 4. The n=5
adapter check verified positive mappings 6 and 10, rejected boundaries 5 and
11, and confirmed all new output/lock paths absent before submission. Review
the complete n=10 prefix before further expansion.

The later grouped-CV nuisance-refit pilot is complete (20 conditional draws
per rho at rho=0/1, plus identities) and does not establish sampling coverage.
The new independent stage is specified in `INDEPENDENT_INFERENCE_PILOT.md`.
Its separate runner is `run_independent_inference_pilot.R` / `.sh`; the old
seed-1--10 calibration runner below remains unchanged.

The new manifest reserves task IDs 1--100, simulation IDs 10001--10100 and
weight-bootstrap B=1000. It retains the original analytic CI and adds a separate
six-row diagnostic inference table. Input/CI/failure/publication/scheduler tests
passed 122 assertions; single and matching-array stubbed launches passed and
the mismatched array was rejected. The current v19 software gates still match.

The initial launch submitted only task 1 as job `18279345_1`, with 40 CPUs,
32 GB and a 24-hour limit; it was first verified PENDING at scheduler priority.
At that point no later task had been submitted. The complete output and six-rho numerical gate must pass before
expansion. The output root is
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
Its `runner_validation.txt` records all gate/fingerprint evidence. The workflow
fingerprint is `9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4`,
and the manifest SHA256 is
`1dc5b92cf549b417a71eadf1cfe4be50bf8a8219e147e171d53f073ff763a16d`.

The first task subsequently entered RUNNING. Two new postprocessing tools are
ready outside the frozen nine-file execution workflow:

- `audit_independent_inference_pilot.R BUNDLE_DIR PACKAGE_LIBRARY OUTPUT_DIR`
  validates exact task/manifest/data/fit structures and 12 numerical checks per
  rho, including arm contrasts, site variance, fold aggregation and signed KKT.
  Its gate binds the exact raw checksum manifest and records dependencies.
- `summarize_independent_inference_pilot.R ROOT N OUTPUT_DIR` accepts only the
  predeclared prefix checkpoints. It requires each seed audit under
  `ROOT/audits/seed_NNNNNN`, verifies its binding to current raw bytes, and
  reports paired MSE and CI-coverage differences with Monte Carlo uncertainty.

The seed auditor passed 78 targeted assertions, and the summary passed 36;
an existing saved fit/data pair also passed the structural compatibility check.
These component tests do not certify a real independent-seed output. Exact
fingerprints are recorded in `postprocessing_validation.txt` under the pilot
root. Missing/unknown hard nuisance counters fail closed, while incomplete
warning capture is explicitly reported as `warning_count=NA`, not zero.
The runner does not provide a complete cross-worker structured warning list;
its Slurm logs remain necessary context for warning review.

### First real seed: execution gate, not statistical confirmation

The raw bundle has 84 rows (six rho settings, 14 methods each), six retained
data/fit artifacts and a six-row independent inference table. All 72 numerical
checks pass; maximum normalized KKT residual is 4.65e-8. Hard nuisance and
weight-bootstrap failure counts are zero. Original estimates and analytic CIs
are unchanged.

| Rho | TATE estimate | Analytic SE | Diagnostic weight-relearning SE |
| ---: | ---: | ---: | ---: |
| 0 | 0.17819069 | 0.02224995 | 0.02345105 |
| 0.5 | 0.18946486 | 0.02294766 | 0.02839791 |
| 1 | 0.17916527 | 0.02507199 | 0.03189765 |
| 1.5 | 0.16560020 | 0.02650100 | 0.03028474 |
| 2 | 0.16560020 | 0.02650100 | 0.02869561 |
| 2.5 | 0.16560020 | 0.02650100 | 0.02751551 |

Truth is 0.2062655638596 for this fixed superpopulation design. Both intervals
contain truth in this one seed at every rho. This is **one** observation per
rho, not six independent coverage replications; each exact-binomial interval
for coverage is [0.025,1]. Empirical SD and Monte Carlo SEs remain unavailable
at n=1. Neither coverage improvement nor RMSE superiority is established.

The repeated log warning is attributable to comparison weights used by
`federated_dr`/`pooled_dr` and their TATE rows, not the primary RoCE correction.
Their persisted diagnostics identify 23 lower clips among 2,000 source weights
(1.15% overall, 2.3% maximum within a source), with limits [0.1,10] and pre-clip
range approximately [0.037999,7.360629]. No upper clips occur. The specific
source label and exact number of fitting calls cannot be inferred from the
saved aggregate diagnostics or repeated log text. All affected comparison
estimates/SEs are finite; retain these diagnostics in subsequent checkpoints.

### Bounded next-prefix array

`run_independent_inference_pilot_array.sh MANIFEST OUTPUT_ROOT PREVIOUS_CHECKPOINT`
is a scheduling-only adapter outside the scientific workflow hash. It verifies
the previous checksum-protected summary, manifest and summary-script identity,
allows only IDs in the next predeclared prefix batch, derives the exact seed
path, and rejects existing output/locks. The n=1 gate permits only tasks 2--5;
it does not authorize all 100 tasks or promote a diagnostic CI.

Two stubbed positive mappings and five rejected task IDs passed without any
fitting. Adapter SHA256 is
`9f84877a8310446e30383abc187b356596c92fd099cd627333e6766d2f72b383`.
The scientific nine-file workflow remains
`9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4`.
The submitted next array is `18290024_[2-5%2]`, initially verified PENDING with
throttle 2 and 40 CPUs per task. Review the complete n=5 checkpoint before any
further expansion; do not overwrite the n=1 audit or summary.

### Retained seed-10003 warning review

The comparison DR diagnostics show 11 lower clips among 2,000 source weights
(0.55% overall, 1% maximum at one source). Log messages identify source-level
counts 1/1000 and 10/1000, but saved aggregate diagnostics do not identify
which count belongs to s1 or s2. These weights serve `federated_dr`/`pooled_dr`
and their TATE rows, not the primary RoCE correction. All related estimates
and SEs remain finite.

Two displayed glmnet warnings concern a binomial class with fewer than eight
observations. They are not bound to exact calls in the saved artifact. A
read-only cell-count review identifies rho=2.5, source s1, A=1/Y=0 as the
strongest candidate: 17 observations in the full source and 12--15 in the
outer training complements, so internal CV training subsets can be smaller.
Target and source-2 cells are much larger. This is a candidate attribution,
not an exact replay or an assertion that the two messages identify two unique
fits. Source/rho/outcome-CV attribution cannot be certified from log order.

All primary degenerate/nonconvergence counters, inference truncation and safety
clip counters are zero, and all six weight-bootstrap failure counts are zero.
The warning review found no implementation blocker requiring suspension of
the existing tasks; warnings and limited capture remain part of the staged
evidence. No model was refitted and no clipping, lambda or simulation setting
was changed to suppress these observations.

### Retained seed-10002 warning review

Seed 10002's six primary TATE rows and bootstrap SEs are finite, its bootstrap
failure counts are zero, and the canonical 72-check audit passes. Comparison
DR diagnostics show 17 lower clips out of 2,000 source weights (0.85% overall,
1.7% maximum within a source), with pre-clip range approximately
[0.04802847,7.481854]. The warning belongs to the already identified comparison
weight helper; repeated log messages are not counted as independent fits.

A lognet small-class warning is also retained. At rho=2.5, source s1 has
17 A=1/Y=0 observations; the five outer-training counts are 14,12,14,15,13.
Target and source-2 minimum outer-training cells are 102 and 87 respectively.
Thus the high-rho source-1 outcome CV is a plausible candidate, but the saved
artifact does not bind the warning to a particular call or internal CV fold.
These counts use source positions resolved through the fit's named weights,
with nonempty observation IDs validated before forming complements. An initial
temporary probe that assumed named outer source-index lists was discarded;
the corrected counts above do not alter any fit, result or numerical audit.
No clipping threshold, lambda rule or data generation setting was changed.

### Seed 10004: a retained noncoverage example, not an optimizer failure

The original analytic CI does not contain truth at rho=0.5,1,1.5 in this seed.
The diagnostic relearned-weight CI contains truth at rho=1.5 but still does
not at rho=0.5,1. These are outcomes from one seed at correlated settings, not
estimates of population coverage and not a reason to change the prespecified
method. All numerical gates, nuisance convergence/degeneracy/line-search and
bootstrap failure checks pass; primary influence tails are not unusually large.

Target-only is 0.2478865 versus truth 0.2062656, an error of +0.0416209 in this
draw. Its analytic SE is 0.0298327 and its CI contains truth. The exact
fold-weighted source increments relative to target-only are:

| Rho | Source 1 increment | Source 2 increment | Net borrowing increment |
| ---: | ---: | ---: | ---: |
| 0 | -0.007562 | -0.011934 | -0.019496 |
| 0.5 | +0.017186 | -0.011249 | +0.005937 |
| 1 | +0.029207 | -0.011951 | +0.017256 |
| 1.5 | +0.025613 | -0.013639 | +0.011974 |
| 2 | +0.010883 | -0.015069 | -0.004186 |
| 2.5 | +0.001955 | -0.015271 | -0.013316 |

The decomposition reconstructs the estimate-minus-target difference to at
most 2.1e-17 using actual fold weights/fractions, not mean weights applied to
pooled estimates. The changed source adds positive error at intermediate rho;
the unchanged source offsets some of it. Source-1 mean fold weights decline
from about 0.364 to 0.005 as rho increases, with corresponding Wald/penalty
increases. Thus the soft filter behaves as implemented but retains borrowing
at intermediate mismatch. A larger diagnostic SE does not remove this point
estimate shift. No cutoff, clipping rule or interval was tuned to make this
seed cover.

Comparison DR clipping is 13/2000 (0.65% overall, maximum source fraction 1.3%),
all lower clips, with finite comparison outputs. Small-class warnings lack
call-level binding; rho=2.5 source-1 A=1/Y=0 has 15 full-site observations and
at least 11 in each outer training complement, so smaller internal CV training
cells are a plausible source. Neither precise warning attribution nor a
complete cross-worker warning count is asserted. The read-only review found
no implementation blocker requiring suspension of task 5.

## API and earlier calibration history

`run_single_simulation()` and `run_simulation_study()` expose the opt-in
argument `n_weight_bootstrap`. Its default is `0L`, which does not call the
bootstrap, consume RNG draws, or alter existing estimates, analytic standard
errors, confidence intervals, or coverage indicators. Positive values must be
at least two.

For the primary one-round soft-penalty TATE estimator, a positive value calls
`estimate_tate_weight_bootstrap()` with weight relearning enabled, no saved
draws, and one bootstrap core. The simulation row retains the analytic standard
error in `se`; bootstrap uncertainty is recorded separately in:

- `variance_weight_relearn_bootstrap`
- `se_weight_relearn_bootstrap`
- `se_fixed_weight_bootstrap`
- `weight_uncertainty_ratio`
- `weight_relearn_n_bootstrap`
- `weight_bootstrap_multiplier`
- `weight_bootstrap_seed`
- `weight_bootstrap_failures`
- `weight_bootstrap_relearn_weights`
- `weight_bootstrap_screening_rule`

The seed is derived without reading or changing `.Random.seed`:

```
((sim_id * 104729 + stream_offset) %% (.Machine$integer.max - 1)) + 1
```

The soft-penalty offset is 13007 and the hard-threshold offset is 17011. `rho`
is intentionally omitted so a same-seed rho group uses common multiplier draws,
matching the paired simulation design. The hard-threshold estimator does not
receive a formal bootstrap SE in this integration; its bootstrap can be run as
a separate post-fit sensitivity diagnostic.

This is a fixed-nuisance, observation-level multiplier bootstrap. Within each
draw, one positive multiplier is generated per original observation and reused
whenever that observation enters an overlapping inner-training fold, its outer
evaluation fold, or either treatment arm. The bootstrap therefore relearns all
fold-specific aggregation weights and captures their dependence on the shared
data, while preserving treated-minus-control covariance. It does not refit the
high-dimensional nuisance models. Its interpretation consequently relies on
the same orthogonality and nuisance-rate conditions as the analytic DML
expansion, and it must be calibrated against a bounded full-refit bootstrap
before being promoted from a diagnostic to the primary standard error.
The implementation supports a fixed numeric aggregation multiplier. Fits
using inner-CV selection of that multiplier are rejected because the CV
pooling path does not yet propagate observation multipliers. Two-round
comparison rows retain missing bootstrap diagnostics in this integration.

Hard thresholding remains a nonregular pretest diagnostic. Under a compatible
source null, its Wald statistic can have a nondegenerate limiting distribution
at the fixed cutoff, so an ordinary n-out-of-n bootstrap need not provide
uniformly valid inference. Hard-threshold bootstrap inclusion frequencies are
useful sensitivity diagnostics, but are not reported as formal confidence
intervals by the simulation driver.

`diagnose_simulation_results()` groups settings by `n_weight_bootstrap` and
fails closed when a requested soft direct-TATE row lacks complete bootstrap
metadata, when variance/SE or relearned/fixed-SE ratio identities fail, when
draw counts or seeds are invalid, or when bootstrap fields appear despite the
diagnostic being disabled. Nonzero bootstrap failures are reported separately.

## Production-dimension calibration pilot

The standalone `run_weight_bootstrap_calibration.sh` / `.R` pair in this
directory runs one seed and all six rho values per Slurm job. It uses the
canonical manifests and locks C1, K=2, p=100, 1,000 observations per site,
five folds, 100 nuisance lambdas, cutoff 1, and 5,000 comparison-bootstrap
draws. Weight-bootstrap B is explicit (200--5,000); the first pilot uses 500.
Only seeds 1--10 are accepted in this diagnostic workflow. Seed 1 is reviewed
before submitting additional seeds. These seeds were already examined in the
v15 experiment, so this is calibration, not an independent confirmation set.

The wrapper requires `ROCE_PROJECT_LIB` and `ROCE_PACKAGE_CHECK_GATE`. Both
installed-package tests and R CMD check must have passed for the current
source fingerprint. The runner also checks the exact manifest mapping,
14-method result schema at every rho, observation-ID contracts, and equality
of the saved fitted objects' point estimates/SEs to the CSV rows.

Successful output is one atomically committed directory containing:

- `results.csv`: all six settings, analytic SE and separate bootstrap fields;
- `diagnostic_qc.csv`: implementation failures and statistical review flags;
- `artifacts.rds`: six rho fits, raw inner/outer components, generated site
  data, task, arguments, and resource plan for subsequent resampling;
- `metadata.txt`: package/workflow/manifest fingerprints and run identifiers;
- `sha256.txt`: payload hashes that can be checked with `sha256sum -c` from
  the output directory.

The first v17 pilot, job `18138523`, failed its artifact-completeness check
because vector and scalar rho formatting produced different keys. It did not
publish results. The group helper and runner now consistently use scalar rho
keys. The fix passed regression tests and a real six-rho smoke at 1,000
observations per site, five folds, 100 nuisance lambdas, and p=4; that smoke
checks the handoff, not production-dimension statistical performance.

The replacement v18 pilot is job `18150037`, dependent on successful jobs
`18149939` (installed tests) and `18149998` (R CMD check). Its output directory is
`results/direct_tate_mc500_b5000/weight_bootstrap_calibration_v18/seed_000001`.
The launch uses 40 CPUs, 32 GB, and at most four positive-rho workers with
five nuisance-CV threads per treatment arm. Runtime and artifact size still
require measurement on this first production-dimension run.

Before statistical expansion, inspect fixed-bootstrap/analytic SE agreement,
weight-relearning SE, source weights, Wald statistics, optimizer diagnostics,
and clipping for every rho. A single successful pilot proves execution and
artifact consistency only. Coverage assessment and a bounded full nuisance-
refit comparison remain necessary before adopting a different primary SE.

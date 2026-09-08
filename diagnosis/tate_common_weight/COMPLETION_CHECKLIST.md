# Full-goal completion checklist

The software-tested outcome CV training-scale correction has been integrated
into workspace source; see `OUTCOME_CV_SCALE_INTEGRATION.md`. Frozen v19
artifacts retain their old scale. This is not a completed inference gate or
authorization to merge results across package versions.

Completion remains unproven. Software checks and the bounded C1/K2 calibration
are prerequisites, not replacements for the full simulation and RHC scope.

Current user constraint: solve the estimator, weight-learning and analytic
inference problem structurally, not with bootstrap CIs or SE inflation.
Historical bootstrap evidence is diagnostic only. The current repair contract
and a verified dual-basis calibration tangent counterexample are recorded in
`STRUCTURAL_TATE_REPAIR.md`. Operational capture integration is paused while
these structural issues are addressed; no new primary method has been adopted.

The bounded training-only target projection selection and its reproducibility
checks are recorded in `TARGET_PROJECTION_VALIDATION_REVIEW.md`. This is a
single-seed, single-outer-fold development check, not a completed inference gate.

## Canonical experiment definition

Use the cutoff-1 manifest universe under
`results/direct_tate_mc500_b5000/production_v35_c1/`. Despite the directory
name, its manifest covers **all C1--C3**, not C1 alone. Direct inspection gives
27,000 primary task keys and 4,500 grouped jobs:

- C1, C2, C3;
- K=2,4,8 source sites;
- rho=0,0.5,1,1.5,2,2.5;
- 500 replications per config/K/rho setting, p=100;
- six rho-specific CSVs per all-or-none config/K/seed commit.

The authoritative cutoff-1 manifest hashes are
`35b8a5d09a59c96d21f56edb2a7fb67ba0dc29e698560766ffca6e479f5d8c43`
(primary) and
`ad4af6e38e2fb1f6ce69a02a9bfdc51ae8897fe6e15174a377e80448a53a8341`
(grouped). The similarly named manifests directly under
`results/direct_tate_mc500_b5000/` also contain 27,000 keys but specify the
older cutoff **2**. They are a different experiment and must not be combined.
K-specific manifests partition the same keys; they add no extra replications.

## Required evidence

| Requirement | Evidence needed | Current status |
| --- | --- | --- |
| Common TATE weights and variance | Exact contrast/covariance/scale tests; documented adaptive-weight conditions; appropriate sampling evidence | Fixed-weight algebra verified; adaptive inference unresolved |
| Reproducible software | Full installed R tests and R CMD check against the final source/package hashes | v18 and grouped-CV v19 passed; repeat if final package changes |
| Warning provenance | Retained warnings with task/method/worker context, explicit capture completeness, and no silent conversion of missing counts to zero | v19 independent-pilot warnings remain in Slurm logs; complete structured cross-worker capture is missing and must be addressed before final production freeze |
| Conditional refit diagnostic | Prespecified grouped-CV identity and 20-draw sets at seed 1, rho=0/1; complete provenance and intermediate audits | Complete: 42 exact bundles/audits; this is not an independent-sampling coverage gate |
| Production freeze | One declared package, workflow, manifest set, nuisance/variance policy and method schema | Final inference/schema not frozen |
| Main simulations | 4,500 grouped commits, 27,000 exact task keys, 54 settings each with 500 replications, no duplicate or mismatched rows | No final-version production task rows; current v18/v19 work is diagnostic only |
| Staged gates | For each of nine config/K blocks: cumulative n=1,5,10,25,50,100,200,300,400,500; implementation and statistical review covering all six rhos | Older-version C1/K2 gates are evidence only for those runs, not final completion |
| Cutoff/inference-radius sensitivity | 1,000 C3/K4 endpoint-rho keys with the prespecified reused grid sidecars, or the explicitly chosen standalone fallback | Final-version outputs missing |
| Fitting-radius sensitivity | C3/K4, rho=0/2.5, M_fit=4/6, M_inf=5, 500 repetitions: 2,000 separate tasks | Final-version outputs missing |
| RHC primary and sensitivities | Matched final package/data/workflow; min numerical audit, conditional complete-1SE fallback, support/pseudo-value/cutoff/reproducibility audits | Local data/cohort verified; final-version analysis missing |
| Final reporting | Re-audited aggregates and C1--C3 TATE figures plus final RHC tables/figures; honest uncertainty and provenance | Not complete |

The 1,000 reused sensitivity sidecars each contain five cutoffs
{1,1.5,2,2.5,3}, plus inference radii {4,6,Inf} at the primary cutoff, for
eight unique grid points. The standalone cutoff manifest is a fallback, not
another 1,000 observations to count on top of those same keys. Fitting-radius
changes require nuisance refits and are separate from inference-radius reuse.

The production method schema must be declared before expansion: older grouped
documentation describes 13 rows while the current bounded calibration includes
14 with the hard-threshold diagnostic. Hard screening is not promoted to the
primary estimator by appearing in a diagnostic CSV. A green gate for a different
schema cannot certify the final output.

The warning-provenance gap is distinct from a numerical failure. In the v19
independent pilot, seeds 10008/10009 passed the exact numerical and payload
audits, but retained comparison DR clipping warnings are not fully captured as
structured conditions. Inspection of `run_direct_tate_rho_group_task.R` and
the rho-group PSOCK path in `R/simulation.R` found no complete structured
warning return at the task/worker boundary. Logging only at the parent process
is not enough for warnings emitted in workers; repeated method/rho messages
must not be mislabeled as independent nuisance fits. Preserve the current
pilot's explicit `warning_capture_complete=FALSE` and archived logs. Address
capture and test numerical/RNG invariance in the final-production version;
do not change the frozen v19 package/workflow during this pilot or relabel its
existing artifacts as having complete warning counts.

The read-only boundary review identifies positive-rho PSOCK workers and nested
treatment-arm/source workers as separate capture points. A future opt-in
condition envelope must return both the unchanged scientific value and raw
warning events from each enabled worker boundary; parent-only handlers are
insufficient. Retain scientific coordinates and an event sequence, leave
unavailable fold/stage context explicitly missing, and distinguish repeated
console propagation from identical messages emitted by different actual fits.
Before adoption, test numerical/RNG invariance across serial, fork and PSOCK
paths, exact nested event counts, warning-free typed output, and preservation
of warnings preceding errors without converting a failed fit to success.

The comparison paths also require explicit worker transport. In
`R/comparison_methods.R`, `.parallel_site_fits` runs independent site AIPW/DR
component fits in workers, and `.resolve_dr_weights_by_site` separately fits
source-to-target density ratios in workers. Capturing only the RoCE
arm/source workers and enclosing rho PSOCK call would still lose comparison
conditions emitted in child processes. Propagate events through the callers
and caches, keeping shared fit events under their original identities rather
than relabeling reused components as new fits. The current production grouped
driver assembles method rows and sensitivity sidecars but has no structured
warning publication/failure payload; that output boundary must also be wired
and tested before declaring complete capture. This inspection does not change
the frozen v19 code or retrospectively supply missing warning provenance.

An isolated source candidate is now available at
`results/direct_tate_mc500_b5000/condition_capture_candidate.1n18nV/`.
Its explicit envelopes now connect crossfit source/arm/refit, comparison-site
and single-simulation/rho-group boundaries. The independent direct-overlay
test total is 202, with zero failures/errors/skips. Bounded serial pipeline
job `18342999` passes both v19-versus-candidate disabled/enabled scientific/RNG
comparisons across three rhos, plus exact warning-message multiset fidelity
(33 original warnings become 33 distinct events, with no escaping warnings).
This uses isolated candidate R functions over tested v19 compiled dependencies,
not an installed candidate. Actual-callsite failure tests and a portable
source/namespace fixture are included, but real parallel fit checks, full
non-CRAN installed tests and production success/failure publication remain
pending at that stage. Standard build/check job `18344636` failed on four fixture-related
assertions in the no-source installed environment; its output is retained.
After repairing only the fixture function lists, all 202 targeted assertions
pass from `/tmp` against the retained installed candidate, without namespace
mutation. Fresh standard check `18345314` now passes with exact `Status: OK`.
Installed-suite job `18346736` also passes: 2,300 repository expectations with
the standard non-CRAN/C2-CV-scale flags, plus 202 capture expectations; 17
additional opt-in repository diagnostics remain skipped and are not counted
as passes. All tracked inputs are unchanged. Enabled real-parallel capture
now passes a bounded installed-package check: six cold-process serial/fork/PSOCK
off/on cases (array 18362596) and audit 18362802 pass 15 checks, with identical
science versus v19 at 1e-12, within-backend on/off RNG equality, 33 matching
typed events and no escaping enabled warnings. This remains p=3, CV-thread=1
evidence, not production-dimensional or statistical validation. Candidate-only
publication/controller helpers additionally pass 101 no-model assertions;
real production-dimensional replay 18372548 subsequently failed closed because
comparison-worker context normalization erased inherited sim_id/rho. No
scientific CSV/commit was published; the original warning payload survives in
the atomic failure bundle. An isolated context-repair source adds six
regression assertions (old implementation fails all six; repair passes).
New R CMD check 18384704, installed suite 18385025 (2,300 repository + 208
capture assertions; 17 other opt-in skips) and topology array 18385031/audit
18385059 (18 checks including inherited coordinates) all pass. Repaired full-
dimensional replay 18385710 completed 0:0 and published six CSVs plus shared
events. Auditor 18385836 failed on automatic CSV type guessing for the all-empty
rho=0 reuse label; a later local audit exposed legitimate NULL warning calls
suppressed by their emitters. Both auditor versions are retained. The final
typed-reader/stage-context auditor passes 22 assertions and real audit_v3:
426 shared scientific fields match exactly, with verified row/event hashes.
Twelve original warning calls remain explicitly unavailable; stage/simulation/
rho context is retained, and warning_capture_complete remains FALSE. This
passes real publication/scientific parity, not full provenance policy or
inference validation. See that directory's
`CANDIDATE_STATUS.md` for hashes and test-scope limitations. It is not deployed;
the warning-provenance requirement remains open. Main-source v19 and its frozen
simulation workflow remain unchanged. The earlier `18340173` crossfit-only
smoke remains separate historical evidence, not an additional MC replicate.

RHC launch software enforcement is now implemented independently of the frozen
pilot: `run_rhc_direct_tate.sh` requires matching installed-test and R CMD check
gates before R or analysis-output creation. Its 27 isolated no-model cases and
an independent rerun pass; the real v19 gates also pass a preflight stopped
before module loading. See `REAL_DATA_FINAL_PREFLIGHT.md`. This closes that
launch-precondition gap, not the pending RHC fitting/inference requirement.

## Diagnostic evidence must not be double-counted

- The v15 common-weight/hard-screen investigation contains 100 C1/K2 commits
  and 600 CSVs. It established important point-estimation and undercoverage
  evidence, but has different package/workflow hashes from v18.
- The old v35 run contains five C1/K2 commits and 30 CSVs under another package.
- The older cutoff-2 run contains ten commits and 60 CSVs under another design.
- v18 weight calibration seeds 1--10 reference the full cutoff-1 manifests but
  execute a separate diagnostic workflow with saved fitted objects. They do not
  turn into formal production rows solely by sharing manifest keys.
- Full nuisance-refit draws are draws conditional on saved datasets, not new
  independent Monte Carlo datasets. Identity draw zero is not a bootstrap draw
  for estimating SD and cannot count toward coverage.

## Communication comparison and reporting boundaries

The inspected canonical manifests contain no current common-weight one-round
versus two-round comparison specification. Existing manuscript text explicitly
labels the older continuous-outcome comparison as archived and descriptive.
That archive must not be relabeled as current common-weight equivalence.
If a current-method communication claim is required in the final manuscript,
its settings, repetitions and gate must be frozen before running that additional
comparison; no new large experiment count is invented here. The absence of that
manifest does not erase any already prespecified primary or sensitivity task.

See `REAL_DATA_FINAL_PREFLIGHT.md` for the real-data execution contract and
`AGGREGATION_OUTPUT_CONTRACT.md` for the distinction between actual foldwise
weights and their reported average. JASA and immutable scientific archives
remain untouched throughout the current diagnostic work.

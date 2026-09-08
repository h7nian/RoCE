# Independent-sampling inference pilot: predeclared protocol

## Purpose and limits

This protocol is fixed before generating or inspecting its simulation results.
The completed grouped-CV pilot used one saved dataset, two rho settings and
20 conditional draws per rho. It established executable, auditable nuisance
refitting; it did not establish confidence-interval coverage. The next pilot
addresses sampling variation across independently generated datasets without
changing the RoCE point estimator, nuisance rule or original analytic CI.

This is separate from the canonical 27,000-key production experiment and all
its sensitivities/RHC requirements. It does not replace any of those tasks.
The final inference policy remains unresolved. A numerically successful pilot
must not be called a passed coverage gate or a proven variance correction.

## Fixed design

One task generates one independent simulation seed and all six paired rho
settings, reusing only the mathematically permitted within-seed rho components.

| Field | Fixed value |
| --- | --- |
| Task IDs | 1--100 |
| Simulation IDs | 10001--10100, task ID + 10000 |
| Configuration | C1, p=100, two source sites |
| Sample size | 1,000 per site, including target; 3,000 total |
| Rho grid | 0, 0.5, 1, 1.5, 2, 2.5; changed source s1 |
| Outcome/estimand | Binary / superpopulation TATE |
| Cross-fitting / nuisance CV | Five outer folds; 100 nuisance penalties; minimum CV loss |
| Fitting / inference bounds | M_tau=M_tau_inference=5 |
| Aggregation | Existing one-round common TATE soft penalty, cutoff 1 |
| Comparison bootstrap | 5,000 draws, unchanged |
| Weight-relearning bootstrap | 1,000 observation-level multiplier draws |
| Methods | Existing six-method input list; existing 14-row output schema including hard-screen diagnostic |
| Software | Tested grouped-CV v19 installation; exact hashes required |

The reserved simulation IDs are disjoint from the canonical production IDs
1--500 and the earlier calibration IDs 1--10. No source-code/workflow references
to these reserved pilot seeds were found in the preflight search. Existing
output directories are never overwritten. This pilot generates new datasets;
its rows are not resampling draws from the saved seed-1 dataset.

The method input list is `one_round_crossfit`, `target_only`, `sample_size`,
`inverse_variance`, `federated_dr`, `pooled_dr`. Hard screening remains a labeled
diagnostic and is not a primary method selected by this protocol.

## Intervals and reported quantities

Original `estimate`, `se`, `ci_lower`, `ci_upper` and `coverage` remain unchanged.
For the primary soft one-round TATE row only, a separate audit table derives

```text
diagnostic_relearn_CI = estimate +/- qnorm(0.975) * se_weight_relearn_bootstrap.
```

Its truth-based coverage is a diagnostic outcome, not a replacement for the
existing coverage field or a claim of valid inference. No resampling bias
correction, favorable-seed removal, cutoff change, variance inflation factor,
or maximum-of-SE selection is permitted in this pilot.

For every rho, report across independent seeds:

- mean error and its Monte Carlo SE, empirical SD, RMSE and paired squared-error
  differences against target-only;
- mean analytic, fixed-weight-bootstrap and relearned-weight-bootstrap SEs,
  their ratios to empirical SD, and interval lengths;
- both prespecified coverage proportions, exact binomial Monte Carlo intervals,
  and paired coverage differences with their Monte Carlo SEs;
- bootstrap failures, nuisance/optimizer diagnostics, warning counts and
  source/fold weights, Wald statistics, discrepancies and influence tails.

Do not pool rho settings as independent observations. Within each seed, rhos,
treatment arms and bootstrap multipliers are deliberately paired. The earlier
full-refit reference used fold-stratified multinomial counts; this API uses
observation-level exponential multipliers. The pilot does not pretend that
those are identical resampling distributions.

## Execution gates

Only task 1 is submitted initially, after local runner/input/derived-CI tests
and the existing v19 installed-test/R CMD check gates pass. An artifact/CSV
identity and complete six-rho numerical audit must pass before expansion.

Prespecified cumulative independent-seed checkpoints are n=1,5,10,25,50,100,
always the corresponding prefix of task IDs. At each checkpoint, distinguish
implementation failures from statistical findings. Retain attempted failures,
warnings and their identifiers; never replace them with fresh seeds. A hard
implementation failure pauses expansion for diagnosis. An unfavorable RMSE
or coverage estimate must be reported with Monte Carlo uncertainty, not hidden
or fixed by post-hoc tuning. Decisions on final inference require a separate
explicit review after the intended evidence is collected.

One seed is an execution gate, not a statistical gate. Small prefixes are
descriptive only. Even n=100 gives substantial Monte Carlo uncertainty for
95% coverage; final production validation and method/theory alignment remain
required. This predeclaration does not authorize silently promoting the
diagnostic interval to the primary paper interval.

## Reproducibility and artifact contract

The new standalone runner explicitly passes `n_weight_bootstrap=1000`; existing
production drivers are not silently altered. It requires an exact manifest,
v19 package/source/software gates, workflow and input fingerprints, and unique
task/output paths. One atomic output bundle per seed contains original rows,
the separate inference audit, numerical QC, saved six-rho data/fits, metadata
and SHA256 checksums. Failed attempts must remain identifiable and cannot
publish a success marker.

Preflight failures before a validated task/claimed slot are identified by the
nonzero Slurm job and its log. Scientific/validation failures publish an atomic
failed-attempt bundle. Publication failures attempt the same failure record
only if the output slot is still absent, and always rethrow the original error;
an existing output is never overwritten. A filesystem failure can prevent even
the failure record from being written, in which case the nonzero job/log remains
the evidence. No such attempt is counted as a successful seed.

Planned output root:
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
The first job uses 40 CPUs and five nuisance-CV threads with the already tested
resource-topology rules. Further concurrency is a scheduling choice made only
after the first execution gate; it must not change statistical inputs.

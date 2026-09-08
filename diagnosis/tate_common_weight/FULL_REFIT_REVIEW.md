# Bounded nuisance-refit diagnostic

Current status: both origin-grouped v19 identity jobs (`18210893`, `18210910`)
passed numerical audits. The first actual grouped rho=0 draw (`18212040`) also
passed row-map, actual origin-partition, moment, variance and optimizer checks;
the first actual grouped rho=1 draw (`18213685`) passed those reviews as well.
The predeclared 20 positive draws per rho are complete and audited. Arrays
`18217670` and `18218871` finished without replacement draws; the four-active-job
ceiling was retained, with rho=1's throttle raised after rho=0 finished. See
`GROUPED_CV_IMPLEMENTATION.md` for current evidence. The historical row-CV
protocol below is retained for interpretation: its draws 2--4 completed, while
its unstarted draws 5--20 were subsequently cancelled after all were verified
unstarted. No output was deleted. Historical row-CV results are not an accepted
inference reference and are never pooled with the corrected protocol.

This workflow compares nuisance-refitted TATE draws with the existing
fixed-nuisance weight diagnostic. It does not change the RoCE estimator,
analytic standard error, primary confidence interval, or manuscript.

## Resampling contract

The input is a checksum-verified frozen v18 FACE C1/K2 calibration bundle.
Each draw samples whole observations with replacement within each site and
original outer fold. It preserves fold sizes, but not treatment-arm counts.
One row map is shared by both treatment arms, all inner/outer training uses,
and all fields. The draw seed does not depend on rho. Duplicate observations
cannot cross original outer/inner cross-fitting folds. Fresh fold views point
to resampled data, and target caches are recreated for every fit.

Only `run_tate_crossfit()` is rerun; unrelated comparison estimators and their
B=5,000 bootstrap are not computed. The same refitted nuisance results support
two paired estimates: newly learned common TATE weights and original foldwise
common TATE weights. Their difference therefore isolates the effect of
relearning aggregation weights within this resampling distribution.

For a matched fixed-nuisance reference, the **same** multinomial row counts
also reweight the saved inner moments and outer pseudo-values through the
existing v18 weight-bootstrap helpers. Each draw thus records four paired
estimates: fixed/refitted nuisances crossed with original/relearned weights.
This avoids changing the resampling distribution when comparing nuisance
refitting. It remains a diagnostic, not a validity proof.

Draw zero retains exactly the original data. It must reproduce the reference
estimate, analytic SE, and fold weights to 1e-10, including reconstruction with
original weights. Positive draws are not interpretable until this identity
gate has passed. Any failure is retained with its message; there is no silent
redraw-until-success or deletion of failed draws.

## Important limitations

- This is conditional on the saved design matrices and original partition,
  not a raw-data full-pipeline bootstrap. It must not be used to claim that
  preprocessing uncertainty in RHC or other DGPs has been included.
- Nuisance-model internal CV still splits row positions, not original IDs.
  Duplicate copies can share a nuisance-CV training/validation pair, although
  they never cross an outer/inner cross-fitting train/evaluation boundary.
  This reruns the existing algorithm on empirical samples but is **not** a
  grouped-origin CV bootstrap; tuning loss can be optimistic. Eliminating this
  requires explicit group-ID plumbing or a separate weighted fitting path.
- The existing fixed-nuisance diagnostic uses site-level exponential weights;
  this comparator uses fold-stratified multinomial counts. Compare nuisance
  refitting against the matched-count reference saved in the same draw, not
  directly against the earlier exponential-bootstrap SEs.
- Soft active-set changes and hard screening can be nonregular. Neither a
  larger bootstrap SD nor successful full refitting establishes valid coverage.
  Independent sampling replications and a prespecified CI rule remain required.
- Sparse-arm failures are part of the attempted-draw record. Conditioning on
  successful draws without reporting failures would change the diagnostic.

## Execution and evidence

`test_full_refit_resampling.R` covers shared row maps, exact partitions,
fresh data references, RNG preservation, input rejection and source-order
equivariance. With seed 1's six saved v18 fits it also reconstructs original
and target-only estimates. The initial 95 assertions passed; adding matched
multinomial counts, all-one weight recovery, and perturbed-count functional
checks increased this to 131 passing assertions before the first full refit.
These are component checks, not a full-refit or statistical validity gate.

`run_full_refit_draw.sh BUNDLE_DIR RHO DRAW_ID OUTPUT_DIR` runs one draw per
Slurm job using five nuisance-CV CPUs and sequential source/arm scheduling.
The shell verifies installed-test and R CMD check gates. The R driver checks
input payload hashes, locked v18 provenance, and package/workflow hashes again
before atomic publication. An existing destination is never overwritten.
Each successful or failed fit attempt records a CSV and an RDS containing
draw details, warnings, row maps (successful fits), and source provenance.

Initial execution gate: seed 1, rho=0, draw=0 (identity). No production
bootstrap CI or full 500-replication expansion is authorized by this gate.

Submitted identity job: `18183278`, five CPUs, 16 GB, four-hour limit;
initially observed PENDING after submission. Output destination:
`results/direct_tate_mc500_b5000/full_refit_calibration_v18/seed_000001_rho_0_draw_0000`.
Submission is not a passed identity gate; the terminal job status, output hashes,
and four paired estimates must still be audited.

The job subsequently completed 0:0 in 02:01:22. Independent auditing verified
both output hashes and all four source-bundle payload hashes. All four paired
estimates equaled the original 0.18914932744483715 exactly. Full-refit SE
0.021690477629616, variance, target-only estimate/SE, fold weights, average
weights and Wald statistics also matched exactly; inclusion matrices were
identical. Matched fixed-nuisance relearned weights differed by at most
5.55e-17. Site row maps were identity vectors 1:1000, multiplicities were all
one, and full outer/inner observation-ID validation passed for both fits.
There were no warnings or failures. The numerical identity gate is therefore
passed; this does not validate the distribution of nonidentity draws.

The point/SE summary is stored at
`results/direct_tate_mc500_b5000/full_refit_calibration_v18_summary_seed1_rho0_identity_v1/`.
It correctly records zero nonidentity draws and no validated inference.
Following the passed identity gate, two bounded jobs were submitted with the
same frozen scripts and five CPUs each:

- `18199028`: seed 1, rho 0, draw 1 (first actual resample).
- `18199029`: seed 1, rho 1, draw 0 (identity for the rho-reuse reference).

Their destinations use the same root and corresponding seed/rho/draw names.
They are submitted work, not completed resampling evidence.

## Positive-draw audit and reversible CV investigation pause

Final queue disposition: after the new v19 grouped-CV software gates and rho=0
production-dimensional identity audit passed, the obsolete row-CV tasks
`18206784_5`--`18206784_20` were individually rechecked as
`PENDING (JobHeldUser)` and cancelled. Slurm confirms all sixteen
`CANCELLED` with elapsed 00:00:00. They were never started; no completed data,
code or audit output was deleted. The exact row-CV protocol archive remains
available for any deliberately scoped historical reproduction. New-protocol
draw counts and the full experiment objective are not reduced by retiring
these old queue entries; use fresh grouped-CV jobs after its own gates pass.

Actual rho=0 draw 1 completed as job `18199028` in 00:46:24 and passed independent
row-map, arm alignment, raw-moment, final variance, Wald/penalty and optimizer
checks. The maximum KKT residual was 9.08e-7. Its plug-in SE increase was traced
to nuisance-prediction changes, not an arithmetic aggregation discrepancy.
The targeted nuisance-CV tests then found actual cross-validation overlap of
duplicate origins (460/800 PS rows and 156/292 OR rows involved). A controlled
grouped-origin sensitivity removes this overlap and changes the identified
extreme AIPW value from -63.95 to -2.70, while also changing fold composition.
This is evidence of an important resampling-CV semantics problem, not proof
that grouping alone establishes valid bootstrap inference.

After the initial computational checks, rho=0 remaining IDs 2--20 were submitted
as array `18206784` with throttle 3 through `run_full_refit_array.sh` (launcher
hash `4203c22dd989fe7a71e0d3b80c55c4772abd26c724617a7bc27e1c64aa0c98a5`).
The scheduling adapter calls the unchanged frozen single-draw implementation;
it verifies both prerequisite audit bundles and does not change its numerical
workflow hash. Following the confirmed CV overlap, exactly pending IDs **5--20**
were reversibly held with `scontrol hold`; all sixteen were verified
`PENDING (JobHeldUser)`. Running IDs 2--4 were left untouched, as were all data
and completed outputs. Do not release the held jobs blindly: first resolve
and validate the resampling-CV protocol and record the disposition of this
original protocol. This pause is triggered by a validation issue, not selection
of a favorable variance ratio. The prespecified 20-draw set is not complete.

The rho=1 identity job `18199029` separately completed 0:0 in 01:26:22. Its
original/refit estimate 0.18808819311560929, SE 0.02539934368091875, variance,
all fold weights and ten Wald values matched exactly. Row maps were identity,
all moment and topology checks passed, and its fresh audit bundle hashes passed.
No rho=1 positive-draw expansion was submitted during the CV investigation.

## Predeclared next resampling checkpoint

Before inspecting either running job's numerical result, the next bounded
checkpoint is fixed as follows. For seed 1 at rho 0 and rho 1, first require
the corresponding successful identity gate and one completed, validated actual
draw. Then use nonidentity IDs **1--20**, retaining all attempts, under the same
frozen package, workflow and deterministic row-map seed. Keep at most four
full-refit draws concurrently active in this initial batch. Do not stop the
draw sequence when a desired variance ratio appears, replace failed IDs with
new seeds, or count draw zero in SD calculations. Failures pause expansion for
diagnosis and remain explicit records.

This 20-draw checkpoint is a computational and distributional pilot, not an
equivalence or coverage gate. Summarize all four paired estimate distributions,
both covariance decompositions, optimizer/nuisance failures, and weight changes.
Small-B uncertainty must remain visible. Any later extension of the draw count
or dataset/rho set must be recorded before that extension is inspected. The
primary variance policy still requires independent-sampling calibration at
adequate Monte Carlo size; even a successful 20-draw nuisance-refit reference
does not settle it.

`summarize_full_refit_draws.R` aggregates an explicit draw-ID set including
identity draw zero; zero is excluded from SD calculations. It verifies payload
hashes, CSV/RDS agreement, source-bundle provenance and point/SE identity, and
retains failed draws. `diagnostic_draw_set_complete` is an execution flag only;
`inference_validated` remains FALSE. The summary's identity scope does not claim
to independently recompute the original fold-weight comparison made by the
frozen draw runner. The standalone synthetic test script covers missing,
duplicate, failed, provenance-mixed and valid-rehashed CSV/RDS-inconsistent
draws. The actual identity bundle subsequently passed this summary audit too.

The summary subsequently gained two explicitly separate covariance scopes.
Writing A for fixed nuisances/original weights, B for fixed nuisances/relearned
weights, C for refitted nuisances/original weights, and D for both refitted:

```text
weight = B - A
nuisance = C - A
interaction = D - C - B + A
D = A + weight + nuisance + interaction
```

`change_D_minus_A` reconstructs Var(D-A). `full_D` reconstructs Var(D), including
Var(A) and **all six** pairwise covariance terms among A and the three changes.
The variance of the change alone is not the increase in total variance; in
particular, baseline-change covariance must not be omitted. Nonconstant
baseline/changes regression tests verify both reconstructions to 1e-12.
With zero or one successful nonidentity draw, decomposition values are NA,
not zero. No existing identity output was overwritten by this summary update.

Frozen diagnostic source hashes for this submission:

- `full_refit_resampling.R`: `895e5e6e0ecc837f41b4dd102e36421eb41e715eb809d34b49e6c58131ee2480`
- `run_full_refit_draw.R`: `a44c0e01b52ecbc440116a9bfa8836b07341f5c373673b7253c1ffd76be0c875`
- `run_full_refit_draw.sh`: `8abd13f243093c23210e2ba4bf24dd00e9b0c957138f3521fffda83ccf3323ea`

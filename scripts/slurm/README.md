# MSI TATE simulation workflow

The production manifest contains only the manuscript's `p = 100` negative-
transfer experiment. It has 27,000 canonical setting-replication rows:

```text
3 configurations × 3 source counts × 6 rho values × 500 replications
= 27,000 setting-replication result keys
```

There are no `p = 50` rows in `manifest_main.csv`. These rows are result keys,
not separately submitted tasks: production does not submit 27,000 Slurm jobs. The
grouped manifest combines the six same-seed rho values for each
configuration/source-count pair and therefore contains exactly

```text
3 configurations × 3 source counts × 500 replications = 4,500 Slurm jobs.
```

Each grouped job independently fits rho zero, then reuses only components
that are mathematically and bitwise invariant across rho. It refits the
changed source's treated-outcome nuisance functions and recomputes the full
TATE aggregation for every positive rho. Strict validators reject any
change in the target site, covariates, treatments, control outcomes, or an
undeclared source. The six outputs retain their original 27,000-row task IDs,
so aggregation and per-setting diagnostics remain unchanged. The older
K-specific manifests are exact partitions of the master and remain available
as a conservative independent-fit fallback; they are not additional
experiments and must not be added to the grouped count.

For K=2 and K=4, grouped production uses five OpenMP threads for each
five-fold nuisance CV path and requests `2K × 5` CPUs. K=8 defaults to two CV
threads and 32 CPUs to improve MSI queueability. The rho-zero reference fits
the two treatment arms concurrently; positive-rho updates refit only the
declared treated source. Compiled kernels default to one CV thread outside the
audited wrappers, preventing accidental nested oversubscription.

The 4,500 grouped rows are still not submitted at once.
`submit_rho_group_direct_tate.sh` defaults to exactly one job with concurrency
one, caps each call at five jobs and concurrency two, and never crosses a
configuration/K boundary. It defaults to MSI's `msismall` partition, while
`ROCE_SLURM_PARTITION` provides an explicit audited override. Inspect each
completed batch before advancing.
Every grouped commit is all-or-none across six primary CSVs, required
sensitivity sidecars, and a final sentinel; a failed commit rolls back files
created by that attempt instead of leaving a task looking complete.
For a versioned `.../raw` output root, the grouped submitter defaults its
sensitivity sidecars to the sibling `.../reused_sensitivity_raw` directory, so
two tested package versions cannot silently share sensitivity files.
The submitter automatically stops at cumulative replication checkpoints
`1, 5, 10, 25, 50, 100, 200, 300, 400, 500`. At each checkpoint it queues one
lightweight dependent audit that reviews all six rho settings, writes their
coverage/RMSE/MCSE and nuisance diagnostics, and creates the gate needed to
advance past that checkpoint. Batches between checkpoints require the most
recent gate, so compute cannot silently outrun implementation review without
creating hundreds of redundant cumulative audit jobs.
Every primary manifest row records `n_bootstrap = 5000`; this controls the
multiplier-bootstrap standard errors of SS, IVW, Federated-DR, and Pooled-DR.
RoCE and target-only continue to use their paired influence-function
variance estimators.

Every result row also records `nlambda_init`, and the aggregation and
diagnostic scripts group on it. This prevents runtime-validation grids (for
example, 50 versus 100 nuisance-CV lambda values) from being silently pooled.
Production manifests remain locked to 100. A one-task smoke or the explicitly
bounded `nlambda_validation` manifest (at most 10 paired seeds) may override it
through the builder's fourth command-line argument. The validation manifest
runs only RoCE and its fold-aligned target anchor because the comparison
estimators do not depend on the RoCE nuisance-path length. The value is
forwarded to the target initial fits and every source initial/calibrated fit;
`run_crossfit()` verifies the value returned by every source worker in every
outer fold before it accepts that fold. It does not alter the independently
fit target-only AIPW anchor, whose propensity and outcome regressions retain
the package-standard 100-point `cv.glmnet` path. Every simulation records and
validates separate SHA-256 fingerprints for the tested package,
the simulation workflow files, and the exact manifest. This distinguishes a
driver revision from a package revision and prevents mixed-workflow outputs
from entering a checkpoint or paired comparison. The paired validation also
requires the anchor to be bit-for-bit identical across FACE nuisance grids.

## Required pre-submission checks

On MSI, run the isolated final-source test gate and make R CMD check depend on
it. Use fresh versioned output paths for every candidate:

```bash
ROCE_AUDIT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_vXX \
  sbatch scripts/slurm/run_package_audit_tests.sh
```

This installs from an immutable allow-listed source stage, enables the C(2)
density-ratio CV scaling audit, runs the complete testthat suite, and writes a
source/package-fingerprinted pass gate. `run_r_cmd_check.sh` independently
builds the same allow-listed tree and writes a pass gate only for exact
`Status: OK`. A production batch requires both gates.
After editing an operational R driver, `parse_slurm_drivers.sh` provides a
lightweight syntax-only gate without reinstalling the package.

Treatment-arm scheduling must be compared on the same TATE estimator,
data seed, nuisance grid, and method set. The reusable audit is
`compare_parallel_reproducibility.sh`; it removes only timing, resource, and
installation provenance fields before requiring every statistical and
diagnostic field to be CSV-exact. Production may enable concurrent treatment
arms only after this exact-reproducibility gate and the complete package test
suite pass. Older timing comparisons that fitted different estimands or method
sets are archived diagnostics and are not production gates.

The corrected seed-1 topology audit used `p=100`, C(3), `K=4`, `rho=0`, a
50-point FACE nuisance grid, and the same TATE plus fold-aligned target
method set on both paths. All 105 compared non-timing fields were CSV-exact.
Concurrent arms took 2,954.00 seconds on 40 CPUs; sequential arms took 4,714.23
seconds on 20 CPUs, a 1.60-fold wall-time speedup. This establishes scheduling
equivalence and informs resource allocation; it does not authorize the
50-point grid, which failed the separate predeclared 50-versus-100 adoption
rule. Production therefore retains 100 nuisance lambdas.

After both gates pass, generate the versioned manifests in one lightweight
Slurm job:

```bash
sbatch scripts/slurm/prepare_mc500_manifests.sh
```

The job creates and audits the 27,000-row main manifest, its exact K-specific
partition, the 4,500-row grouped manifest, the one-row smoke manifest, and the separately scheduled cutoff
(1,000 tasks, retained only as a standalone fallback) and fitting-radius
(2,000 tasks) sensitivity manifests. Every
statistical setting contains exactly `sim_id = 1:500`; every main and smoke
row has `p = 100` and `n_bootstrap = 5000`. The sensitivity task counts are
not part of the 27,000-task primary array and are not submitted automatically.
After all structural checks pass, the audit atomically writes
`manifest_audit_passed.txt` with hashes of the master, grouped, and three
K-specific manifests. The grouped production submitter verifies those hashes,
p=100, 500-replicate, and B=5000 declarations before calling Slurm.

The same job then builds the shared-shift family (HISTORY #0009) under
`results/direct_tate_mc500_b5000/shared_shift/` with the same file names: a
9,000-row main manifest (`experiment = shared_shift`, C(1), K = 2/4/8, the six
rho values read as a shared log-odds shift of both arms of `s1`,
`deviation_mechanism = both_arms`), its K-specific partition, and a 1,500-row
grouped manifest, audited with `audit_mc500_manifests.R <root> 500 5000
shared_shift`. Every manifest carries `deviation_mechanism` as its last
column; the grouped submitter, task runners, row annotation, checkpoint audit,
and aggregation read it, and the grouped submitter derives the group count
from the primary manifest instead of assuming 4,500. To submit that family,
pass its root as the manifest root and `.../shared_shift/raw` as the output
root; its result root must hold its own reuse-equivalence gate (below).

The choice of multiplier-bootstrap count can be checked independently of
nuisance fitting:

```bash
sbatch scripts/slurm/benchmark_bootstrap_replicates.sh
```

This compares 1,000, 2,000, and 5,000 draws on a fixed five-site influence
design, writing runtime and Monte Carlo stability summaries under
`bootstrap_benchmark/`. It is a one-task computational benchmark, not a
simulation replicate or a source of paper results.

The 10-run benchmark selected 5,000 draws: the relative Monte Carlo standard
deviation of the bootstrap SE fell from 1.86% at 2,000 draws to 1.16% at 5,000
draws, while mean runtime per five-site bootstrap increased from 0.49 to 1.18
seconds. This cost is negligible relative to the p=100 nuisance fits, so the
more stable value is used in the B5000 production manifests.

For the broader package metadata, namespace, dependency, documentation, and
compiled-code audit, submit:

```bash
sbatch scripts/slurm/run_r_cmd_check.sh
```

This builds the exact working tree and keeps all check artifacts under
`results/direct_tate_mc500_b5000/package_check`.

The array scripts load a pinned project-local installation. Production uses
`Rlib_current` by default; a validation run may set `ROCE_PROJECT_LIB` to an
isolated tested installation such as `Rlib_final_audit`. This prevents a
long-running batch from silently using a different checkout or user-library
version.

Before a production batch, generate and run exactly one `p = 100`, `K = 4`,
C(3), `rho = 0` smoke task:

```bash
Rscript scripts/slurm/build_direct_tate_manifest.R \
  results/direct_tate_mc500_b5000/manifest_smoke_single.csv smoke
ROCE_BATCH_SIZE=1 ROCE_MAX_CONCURRENT=1 \
ROCE_NUISANCE_CV_THREADS=5 \
ROCE_TASK_VERBOSE=1 \
  scripts/slurm/submit_direct_tate.sh \
  results/direct_tate_mc500_b5000/manifest_smoke_single.csv \
  results/direct_tate_mc500_b5000/smoke_final_direct_tate_raw 40
```

The bounded submitter automatically writes this smoke's reused-sensitivity
sidecar to `smoke_final_direct_tate_raw_reused_sensitivity/`, separate from the
production `reused_sensitivity_raw/` directory. This separation prevents the
smoke's `sim_id = 1` from entering the 500-replicate cutoff/truncation
summaries. Earlier smoke outputs remain archived under their original
directories and are never overwritten or admitted to this final gate.

The final `40` requests four source-site workers for each of the two concurrent
treatment-arm fits, with five CV-fold threads per worker. Using `20` remains a
valid lower-resource fallback, but it fits the two arms sequentially and
therefore does not exercise the fully parallel production path.
`ROCE_TASK_VERBOSE=1` is recommended for this single smoke so its log
identifies the active fitting stage; production remains quiet by default to
keep array logs compact.

For a separate runtime-validation smoke with 50 nuisance lambda values, retain
the production smoke and use a distinct manifest and output directory:

```bash
Rscript scripts/slurm/build_direct_tate_manifest.R \
  results/direct_tate_mc500_b5000/manifest_smoke_nlambda50.csv smoke 1 50
ROCE_MANIFEST=results/direct_tate_mc500_b5000/manifest_smoke_nlambda50.csv \
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/smoke_nlambda50_raw \
ROCE_TASK_VERBOSE=1 \
  sbatch --array=1-1%1 scripts/slurm/run_direct_tate_array.sh
ROCE_EXPECT_NLAMBDA=50 \
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/smoke_nlambda50_raw \
  sbatch scripts/slurm/audit_direct_tate_smoke.sh
```

The older multi-row smoke manifests are retained as archived artifacts; they
are not the pre-production smoke command.

For a bounded paired-seed grid study, generate distinct manifests, for example:

```bash
sbatch scripts/slurm/prepare_nlambda_validation.sh
```

Each contains only seeds `1:5` for C(3), `p=100`, `K=4`, and `rho=0`. Submit
these through the same bounded array wrapper with at most two concurrent tasks.
Do not change the production manifest from 100 based on a one-seed smoke; use
the paired study to compare TATE estimates, standard errors, weights,
Wald statistics, convergence diagnostics, and runtime.

The completed 20-versus-50 study is an archived, disqualified experiment: its
largest paired estimate change was 0.571 reference SE and its mean SE ratio was
1.073, both outside the predeclared tolerances. It must not be used to unlock
production. The active decision is the provenance-isolated 50-versus-100 study
described below; production remains at 100 unless all five of those pairs pass.

After both same-seed grid runs are complete, compare their
point estimates, standard errors, TATE diagnostics, and elapsed times
in a lightweight Slurm job:

```bash
sbatch scripts/slurm/compare_nlambda_validation.sh \
  results/direct_tate_mc500_b5000/nlambda20_validation_raw \
  results/direct_tate_mc500_b5000/nlambda50_validation_raw \
  results/direct_tate_mc500_b5000/nlambda20_vs_50_validation
```

This check fails on mismatched seeds, DGP metadata, method sets, estimand
truths, tested-package libraries, or package fingerprints, and on nuisance
nonconvergence, any recorded density-ratio line-search event, a final
density-ratio coefficient at or above 90% of the numerical parameter bound,
or inference safety clipping. This paired-grid adoption rule is deliberately stricter than
the production checkpoint policy: it rejects a computational shortcut after
any such event, even when the selected fit eventually converged. The immutable
library and fingerprint are written into every paired
summary. Candidate lambdas excluded by the bounded CV convergence rule are
counted for review but do not invalidate a run when the selected candidate
converged in every fold. Production manifests remain locked to 100 until this
bounded paired study is complete.

Any faster FACE nuisance grid is adopted only after all five paired seeds
complete and, relative to its reference grid, its mean absolute estimate
difference is at most 0.10 reported SE, its maximum is at most 0.25 SE, its
mean SE ratio is within 0.98--1.02, and its median runtime ratio is at most
0.80. Any final-fit implementation failure disqualifies it. These thresholds
are fixed before reading the paired results; the exercise is a bounded
numerical/runtime equivalence study, not a coverage experiment.

Generate seed 1 with the same bounded wrapper; this is also the recovery path
if a pre-provenance pilot was archived without being admitted to validation:

```bash
scripts/slurm/submit_nlambda_validation_pair.sh \
  results/direct_tate_mc500_b5000 1
```

After the audit for seed 1 passes, advance exactly one paired seed at a time:

```bash
scripts/slurm/submit_nlambda_validation_pair.sh \
  results/direct_tate_mc500_b5000 2
```

The wrapper refuses gaps or existing outputs, submits only the two same-seed
simulation jobs (one per grid), and attaches one lightweight after-success
comparison. A successful comparison atomically writes an implementation-audit
sentinel; the wrapper refuses the next seed until that sentinel exists. Set
`ROCE_SUBMIT_DRY_RUN=1` to inspect the exact pair without submitting it. Do
not advance to the next seed until the statistical differences are also
reviewed.

Before choosing the 50-point candidate, create a fresh, provenance-isolated
50-versus-100 root and reproduce seed 1 on both paths with the same tested
package, workflow, OpenMP thread count, and DGP seed:

```bash
VALIDATION_ROOT=results/direct_tate_mc500_b5000/nlambda50_vs100_v14_openmp
sbatch scripts/slurm/prepare_nlambda_validation.sh \
  "${VALIDATION_ROOT}" 50 100
```

After that preparation job completes, submit the first pair from the same
shell (or set `VALIDATION_ROOT` again):

```bash
ROCE_PROJECT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_20260814_v14_cvcoverage \
ROCE_NUISANCE_CV_THREADS=5 \
  scripts/slurm/submit_nlambda_validation_pair.sh \
  "${VALIDATION_ROOT}" 1 50 100
```

The preparation job creates two five-row manifests, while this first submitter
uses only seed 1. It submits exactly two one-element arrays (one per grid),
requests 40 CPUs for each five-thread task, and attaches a provenance-checked
50-versus-100 comparison. Existing outputs are never overwritten, so later
validation seeds cannot be silently mixed into the seed-1 gate.

If seed 1 satisfies the predeclared numerical thresholds, advance the 50-vs-100
study one paired seed at a time:

```bash
scripts/slurm/submit_nlambda_validation_pair.sh \
  "${VALIDATION_ROOT}" 2 50 100
```

The same audit sentinel and no-gap rules apply through seed 5. The 50-point
grid is eligible for production only after all five same-seed comparisons pass
the mean and maximum estimate-difference, SE-ratio, runtime, convergence, and
provenance criteria. If seed 1 already violates a criterion that cannot recover
(in particular the maximum estimate difference), do not spend four more
paired runs. The 100-point production lock remains in force until this formal
five-seed decision is complete.

After the single smoke task finishes, run its implementation audit:

```bash
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/smoke_final_direct_tate_raw \
sbatch scripts/slurm/audit_direct_tate_smoke.sh
```

This checks the exact method/estimand schema, arithmetic identities, fixed
`p=100`, C(3), `K=4`, `rho=0` metadata, and the RoCE weight, Wald,
optimizer, nuisance-convergence, and clipping diagnostics. One replicate is an
implementation gate only; it is never interpreted as a coverage or RMSE study.
Only after every check passes does the audit atomically write
`direct_tate_smoke_audit_passed.txt`. The production submitter requires this
gate and verifies that its package fingerprint, workflow fingerprint, and
`nlambda_init` match the proposed production batch. It also matches the smoke
audit and production-submitter script hashes, so a stale gate from another
candidate or an earlier operational policy cannot unlock submission. Hashes
of the primary smoke CSV and its reused-sensitivity sidecar also prevent a
previous pass file from authorizing altered or incomplete smoke artifacts.

## Exact grouped-reuse gate and small bounded submission

Before production, run one operational p=100, C(3), K=4, rho=2.5 independent
task and its corresponding six-rho grouped task with the final tested package.
`audit_rho_reuse_equivalence.R` requires every shared non-timing statistical
and diagnostic CSV field to be bit-for-bit identical and writes the reuse gate.
This is an implementation/runtime audit, not Monte Carlo evidence.

The locked setting is selected by `ROCE_REUSE_CONFIG`, `ROCE_REUSE_K`, and
`ROCE_REUSE_RHO` (defaults C3 / 4 / 2.5). The shared-shift family audits its
own both-arm reuse with `ROCE_REUSE_CONFIG=C1` on its manifest root; the gate
records the `deviation_mechanism`, and `submit_rho_group_direct_tate.sh`
requires the gate under the family's result root
(`<result root>/rho_reuse_equivalence_final/rho_reuse_equivalence_passed.txt`)
with the mechanism of the primary manifest.

The final same-package audit used `p = 100`, C(3), `K = 4`, `rho = 2.5`,
seed 1, 100 nuisance penalties, and 5,000 bootstrap draws. All 137 shared
non-timing fields were exact. The six-rho grouped job took 3,666.48 seconds;
using the 2,866.29-second independent endpoint task as the same-package
benchmark gives an estimated six-independent-run/grouped speed ratio of 4.69.
The audit gate is `production_v34_c1/rho_reuse_equivalence_v35/audit/
rho_reuse_equivalence_passed.txt`.

After that gate and both package gates pass, the first production call is:

```bash
ROCE_BATCH_SIZE=1 \
ROCE_MAX_CONCURRENT=1 \
ROCE_PROJECT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_vXX \
ROCE_PACKAGE_CHECK_GATE=results/direct_tate_mc500_b5000/package_check_vXX/r_cmd_check_passed.txt \
scripts/slurm/submit_rho_group_direct_tate.sh \
  results/direct_tate_mc500_b5000 \
  results/direct_tate_mc500_b5000/raw
```

For subsequent batches, omit `ROCE_BATCH_START`: the submitter scans grouped
commit sentinels in order and starts at the first incomplete seed-setting job.
Increase the batch only after reviewing the prior batch. The caps are 500
grouped jobs per call and 50 concurrent array tasks (HISTORY #0010); the
checkpoint ladder truncates every call at the next rung regardless. Set
`ROCE_SETTING=C1:K2` (config and K) to scope a call to one config/K block so
that independent blocks, each walking its own ladder, run concurrently; without
it the submitter walks the grouped manifest in order and requires the previous
block's n=500 gate before starting the next block. The legacy `submit_main_direct_tate.sh` remains
available only as an independent-fit fallback and retains its own small-batch
gates.

The grouped submission wrapper requests 12G, 16G, and 32G for K=2, K=4, and
K=8, respectively. The larger K=8 request matches its 16 simultaneous source/
arm workers while retaining the queue-friendly two-thread CV topology.
Override this only when needed with
`ROCE_MEMORY_PER_TASK`. `ROCE_PROJECT_LIB` must identify an isolated
RoCE installation that passed the audit suite. Both the submitter and array
runner reject a missing/incomplete installation instead of falling back to an
older package on the system library path. Every new task records that library
and a path-independent SHA-256 fingerprint of the content hashes for its
installed DESCRIPTION, shared object, and R lazy-load database. Smoke,
checkpoint, and full-study diagnostics use the same
shared provenance validator: they reject missing, mixed, malformed, or
unexpected package metadata before interpreting coverage or RMSE.
The same readers fail if a production row is not labeled as the p=100 FACE
superpopulation DGP or if its heterogeneity label disagrees with whether
\(\rho\) is zero, preventing scientifically different rows from being pooled.
Submission resolves `Rlib_current` to its immutable versioned target before
exporting the job environment, so a later symlink update cannot change a
queued batch's package. The bounded-array helper also caps the requested row
range at the current experiment/configuration/p/K/rho/cutoff boundary, so even
an explicitly enlarged batch cannot silently mix two statistical settings.
The grouped production submitter rejects more than five jobs per call or more
than two concurrent elements. The independent-fit fallback uses its separate
25-task/five-concurrent safety caps.

Repeat for only one `K` value at a time.  A completed task writes
`task_<global-task-id>.csv` atomically.  Existing task files are never
overwritten; rerunning a completed array index exits successfully after a
skip message. Wait for the current bounded batch and its checkpoint audit
before invoking the submitter again; the missing-output scan is not a queue
reservation mechanism. Set `ROCE_SUBMIT_DRY_RUN=1` to print the resolved
manifest rows, array specification, resources, and immutable package library
without calling `sbatch`.

## Aggregation and per-setting checks

After a batch or the full study completes:

```bash
Rscript scripts/slurm/aggregate_direct_tate.R \
  results/direct_tate_mc500_b5000 500

Rscript scripts/slurm/diagnose_direct_tate_settings.R \
  results/direct_tate_mc500_b5000 500 \
  results/direct_tate_mc500_b5000/manifest_main.csv
```

Diagnostics are computed separately for TATE and retained legacy
potential-outcome-mean rows.  Every configuration/source-count/rho/method
combination is checked for:

- missing, duplicated, non-finite, or invalid-SE replications;
- empirical coverage and its nominal Monte Carlo band;
- bias and its Monte Carlo standard error;
- empirical SD, mean reported SE, and their ratio;
- RMSE and the exact bias/variance moment identity;
- common-random-number paired MSE differences between TATE and
  target-only plus every nonadaptive TATE benchmark, with MCSE, standardized
  difference, and the fraction of replicates favoring TATE;
- nuisance non-convergence counts, maximum final-fit iterations, final
  parameter-update/stopping-threshold ratios, density-ratio line-search
  failures, and maximum initial/calibrated density-ratio coefficient
  magnitudes;
- nuisance-CV candidate fold fits and penalty values excluded because they did
  not converge with finite validation loss in every inner fold, reported as a
  review signal separately from final-fit nonconvergence;
- finite-sample density-ratio support-floor activation counts and the largest
  applied floor, reported separately from optimizer non-convergence;
- minimum site-wide and target-only treatment-arm/outcome cell counts, plus
  the number of cells below eight, so sparse binary designs are distinguished
  from numerical failures;
- for Federated-DR and Pooled-DR, the total and largest site-specific fraction
  of normalized density-ratio weights clipped to the declared bounds, together
  with the pre-clipping minimum and maximum;
- maximum absolute unconstrained source weight, with values above 10 flagged
  for review rather than silently truncated or reset;
- aggregation-optimizer iteration counts and any finite-sample PSD ridge used
  to regularize a noisy plug-in covariance objective;
- the fraction of arm-specific source density logits truncated at
  \(M_{\tau,\mathrm{inf}}=5\), their maximum absolute pre-truncation value,
  and any additional numerical safety clipping;
- normal-reference coverage, which helps distinguish centering bias from
  standard-error miscalibration.

The archived `p = 100` summaries can be audited without rerunning any nuisance
fit:

```bash
sbatch scripts/slurm/diagnose_legacy_p100_summaries.sh
```

This writes immutable legacy setting diagnostics and the `rho = 0` RMSE
ranking under `results/direct_tate_v1/legacy_p100_diagnostics/`; these outputs
remain separate from the TATE production summaries.

The C(3), `rho = 0` centering diagnosis has a separate oracle-nuisance
decomposition that does not refit source models or enter the production array:

```bash
ROCE_REMAINDER_CONFIG=C3 sbatch scripts/slurm/diagnose_target_remainder.sh
```

It runs the exact `p = 100`, `K = 4`, five-fold target estimator for 500
replications and compares fitted/fitted AIPW with three counterfactual
diagnostics: true outcome/fitted propensity, fitted outcome/true propensity,
and both nuisances true. This isolates the finite-sample doubly robust product
remainder from TATE aggregation and variance estimation. Results are written
to a fresh output directory atomically, include the immutable package
fingerprint, and are never overwritten. The provenance-hardened rerun under
`c3_target_remainder_provenance_v4/` reproduces all 15 core columns of the
archived `c3_target_remainder/` raw file byte-for-byte.

The distinct question of why a non-adaptive method can have smaller RMSE at
`rho = 0` has a lightweight DGP-level diagnostic:

```bash
sbatch scripts/slurm/diagnose_rho0_site_tate.sh
```

This performs no nuisance fitting. For `p = 100` and `K = 2, 4, 8`, it
integrates each site's true marginal binary-outcome treated mean, control mean,
and risk difference under the site-specific skew-normal covariate law while
holding the conditional log-odds treatment shift fixed. The output therefore
quantifies both the arm-scale and TATE-scale non-collapsibility/covariate-shift
differences that SS and IVW trade against their variance reduction. It uses
only the four nonzero outcome coordinates,
which is distributionally identical to the `p = 100` outcome signal and avoids
allocating an unnecessary 100-column reference matrix. Results are written
atomically under `rho0_site_tate/` and record the tested package fingerprint.

`flagged_setting_methods.csv` is a review queue, not an automatic declaration
of a code error.  Coverage flags must be interpreted together with bias,
SE/SD, and normal-reference coverage.

Production advances within one setting through initial pilot checkpoints at
1, 5, and 10 completed replications, followed by cumulative checkpoints at
25, 50, 100, 200, 300, 400, and 500. Start with the default one-task batch,
then use batches of at most five jobs (maximum concurrency two) and audit the
exact completed prefix at every gate. Only
after those pilots pass implementation QC and their statistical review queues
have been inspected should a later call use the full hard cap of five jobs and
two concurrent elements. At each listed boundary the grouped submitter
automatically queues `run_rho_group_checkpoint_audit.sh`, which audits all six
rho values; no separate per-rho checkpoint submission is required. The older
`audit_direct_tate_checkpoint.sh` entry point remains available only for the
independent-fit fallback.

The checkpoint driver requires every expected task file, verifies the exact
13-row method schema and task-to-seed mapping, and fails on implementation
problems such as invalid standard errors, CI arithmetic, nuisance
nonconvergence, density-ratio coefficients reaching 90% of their
numerical hard bound, optimizer failures, internally inconsistent cell or
density-ratio clipping metadata, or inconsistent bootstrap metadata. Sparse
cells and valid nonzero DR clipping are retained as statistical-review signals
rather than silently classified as implementation failures.
Rejected intermediate line-search steps are retained in the review queue; a
fit that ultimately reports convergence ended on a successful complete KKT
scan, so that historical counter is not conflated with unresolved
nonconvergence. The shared classification is defined in
`simulation_qc_policy.R` and reused by checkpoint and final-summary drivers. It
writes a separate immutable timestamped audit directory containing the full
TATE coverage/RMSE ranking, paired direct-versus-benchmark MSE audit, and a
statistical review queue. Only after all
implementation checks pass does it atomically create
`implementation_audit_passed.csv`, which records the exact package, workflow,
manifest, and cumulative replication count reviewed. The presence of reports
without this pass gate is not sufficient to advance the array. Coverage outside
its nominal Monte Carlo band, material centering bias, SE/empirical-SD
mismatch, or extreme weights require inspection before the next batch but are
not silently relabeled as code failures. Use cumulative audit counts 1, 5, 10,
25, 50, 100, 200, 300, 400, and 500; only then advance to the next contiguous setting in
the manifest.

The RHC driver uses 5,000 comparison bootstrap draws and the pre-specified
empirical support screen. It requires a passed `audit_tests_passed.txt` in the
chosen library plus an explicit `ROCE_PACKAGE_CHECK_GATE`, with matching
installed-package, source, test-suite and audit-driver fingerprints. Missing,
stale or duplicate required fields stop before R/model execution and analysis
output creation. These software gates do not settle the pending final RoCE
inference policy. After that policy and the final package are frozen, replace
the version placeholders below with their actual matching paths and run the
minimum-loss nuisance rule in a new output directory:

```bash
ROCE_PROJECT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_vXX \
ROCE_PACKAGE_CHECK_GATE=results/direct_tate_mc500_b5000/package_check_vXX/r_cmd_check_passed.txt \
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/rhc_final_vXX_min \
sbatch scripts/slurm/run_rhc_direct_tate.sh
```

This writes to the explicitly chosen output root, excludes both raw
components `Medicare & Medicaid` and `No insurance`, and retains the four raw
strata Private, Medicare, Private & Medicare, and Medicaid without further
collapsing. The existing `results/direct_tate_v1/rhc/` and other `v1` RHC
directories remain immutable archives. They predate the final TATE
support screen and B5000 comparison inference and must not be copied into the
manuscript as the primary result.

If the `lambda.min` nuisance fits fail the strict audit (including a selected
coefficient at 90% or more of the hard numerical bound), run the pre-specified
stronger-regularization analysis in a new directory:

```bash
ROCE_PROJECT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_vXX \
ROCE_PACKAGE_CHECK_GATE=results/direct_tate_mc500_b5000/package_check_vXX/r_cmd_check_passed.txt \
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/rhc_final_vXX_1se \
ROCE_NUISANCE_LAMBDA_RULE=1se \
sbatch scripts/slurm/run_rhc_direct_tate.sh
```

Never overwrite the failed candidate; compare estimates, standard errors,
coefficients, line-search diagnostics, and support diagnostics side by side.
In the historical v35 analysis, the `lambda.min`, `c = 1` candidate is retained
at `production_v35_c1/rhc_primary/` but fails because one control-arm Medicaid
initial density-ratio coefficient reaches the hard bound. The fully propagated
`lambda.1se`, `c = 1` run at `production_v35_c1/rhc_1se/` passes every strict
audit and was that version's primary RHC result, not the pending final-version
analysis. Its pre-specified cutoff grid is
obtained by reaggregation of those same audited nuisances, not by refitting or
choosing a cutoff from the reported confidence intervals.

Historical note: `rhc_supported_1se_v14_openmp_final/` was produced before the
nuisance-rule argument was propagated into the target propensity and outcome
fits.  It therefore means **source 1SE / target min**, not full 1SE.  Its
independent rerun in `rhc_supported_1se_v14_openmp_repro_exact/` reproduces all
739 common scientific fields exactly (RoCE estimate 0.04522971, SE
0.02123078), so the artifact is retained under that precise label.  A result
may be called **full 1SE** only when it is generated by the corrected rule-
propagation implementation and its metadata, target fits, and source fits all
record `nuisance_lambda_rule = 1se`.  Never substitute the historical artifact
for that corrected sensitivity.

The earlier all-site construction collapses the small `No insurance` stratum
into the smallest retained insurance site. The support audit showed that
screening the affected fourth source requires excluding both of its raw
components, which is why the support-screened analysis above is primary. The
all-site analysis remains available only as a separately named sensitivity:

```bash
ROCE_PROJECT_LIB=results/direct_tate_mc500_b5000/Rlib_final_candidate_vXX \
ROCE_PACKAGE_CHECK_GATE=results/direct_tate_mc500_b5000/package_check_vXX/r_cmd_check_passed.txt \
ROCE_OUTPUT_ROOT=results/direct_tate_mc500_b5000/rhc_final_vXX_all_sites \
ROCE_RHC_TOTAL_SITES=5 \
ROCE_RHC_EXCLUDE_SITES='' \
sbatch scripts/slurm/run_rhc_direct_tate.sh
```

The site count and exclusions are parsed by the common RHC runner and recorded
in the result metadata. Every output directory must be new; the runner refuses
to overwrite result artifacts.

Before paying for the full sensitivity run, the focused probe below refits only
the archived control-arm outer-fold/source paths that failed under
`lambda.min`, using the identical fold construction and `lambda.1se` cache:

```bash
ROCE_PROJECT_LIB=results/direct_tate_v1/Rlib_current \
sbatch scripts/slurm/probe_rhc_initial_dr_1se.sh
```

The probe is a numerical stability gate, not a replacement for the full RHC
analysis or its independent audit.

After installing the support-floor implementation, the same focused driver
can verify the formerly failing `lambda.min` paths without duplicating the
fold-construction code:

```bash
ROCE_PROJECT_LIB=results/direct_tate_v1/Rlib_final_audit \
ROCE_OUTPUT_ROOT=results/direct_tate_v1/rhc_supportfloor_probe \
ROCE_NUISANCE_LAMBDA_RULE=min \
ROCE_RHC_PROBE_FILENAME=rhc_initial_dr_supportfloor_probe.csv \
sbatch scripts/slurm/probe_rhc_initial_dr_1se.sh
```

Before copying the fresh result into the manuscript, run the reusable audit:

```bash
sbatch scripts/slurm/audit_rhc_direct_tate.sh
```

The primary audit also requires the exact two-level support exclusion and
5,000-draw bootstrap metadata. It fails fast on pseudo-value point-estimate or
variance inconsistencies, invalid CIs, weights, Wald summaries, cross-arm
covariance identities, nuisance nonconvergence, or non-finite clipping
diagnostics. It writes the machine-readable `rhc_direct_tate_audit.csv` and
comparison audit into a fresh, atomically committed
`audits/<timestamp>_job<id>/` directory beneath the RHC result root. Both
passing and failing audits are retained with a corresponding gate file, so a
later run cannot overwrite the evidence from an earlier run. To audit a
deliberately different sensitivity, set `ROCE_OUTPUT_ROOT`,
`ROCE_EXPECT_RHC_EXCLUDED_SITES`, and `ROCE_EXPECT_N_BOOTSTRAP` explicitly.

After the audit and isolated package tests pass, render the primary RHC forest
plot from the audited method CSV without refitting any nuisance model:

```bash
sbatch scripts/slurm/render_rhc_direct_tate_figure.sh
```

The renderer validates the exact six-method schema, reserves vermillion only
for RoCE, uses the enlarged manuscript font, and atomically installs distinct
PDFs in the result, `docs/figures/`, and `overleaf/figures/` directories. It
does not overwrite the archived arm-wise RHC figure.

After the strict primary audit succeeds, reuse the fitted arm-specific
nuisances for the pre-specified cutoff sensitivity:

```bash
sbatch --dependency=afterok:RHC_AUDIT_JOB_ID \
  scripts/slurm/summarize_rhc_cutoff_sensitivity.sh
```

This post-processing neither refits nuisance models nor reads individual-level
RHC records. The `c = 1` row must reproduce the primary TATE estimate,
standard error, and common weights to numerical tolerance; the script fails
before writing output if that invariant is violated. It also requires the RDS
to match the tested package library and fingerprint, and records that immutable
provenance in both cutoff-sensitivity tables.

To inspect every arm/outer-fold/source/inner-fold nuisance fit without
changing the fitted result, run:

```bash
sbatch scripts/slurm/extract_rhc_fit_diagnostics.sh
```

This writes `rhc_direct_tate_nuisance_detail.csv` plus a pseudo-value schema
and finiteness audit in the same output directory.

When density-ratio coefficients approach the numerical boundary, distinguish
an optimizer defect from empirical support failure with:

```bash
sbatch scripts/slurm/audit_rhc_initial_dr_support.sh
```

The support audit checks each source-arm training design for constant features
whose target--source moment gap exceeds the selected L1 penalty; in that case
the finite-sample initial density-ratio objective has an unbounded
intercept/feature direction.

After all expected setting-method combinations pass implementation QC, submit
the figure job:

```bash
sbatch scripts/slurm/render_direct_tate_figures.sh
```

For a versioned result root, set `ROCE_RESULT_ROOT` to that root and
`ROCE_MANIFEST_ROOT` to the immutable directory containing
`manifest_main.csv`; the figure job does not require copying the manifest into
the result directory.

It aggregates and diagnoses the raw results again, renders only the `p = 100`
TATE panels with the enlarged plotting theme, and installs the three
C(1)--C(3) PDFs under `*_direct_tate.pdf` names in both `docs/figures/` and
`overleaf/figures/`. The archived arm-wise PDFs keep their original filenames
and are also copied to `results/direct_tate_v1/legacy_figures/`; rendering a
TATE figure never overwrites them. The job stops before rendering when
any expected replicate or setting-method combination is missing, or when an
implementation-level diagnostic fails. Coverage, bias, SE/empirical-SD, and
RMSE flags remain statistical review signals: they are written to the audit
tables for setting-by-setting review but do not by themselves suppress a
complete figure.

For a targeted convergence check, `profile_face_fold_components.sh` accepts
`ROCE_NUISANCE_MAX_ITER`. Cross-validation retains its internal iteration
cap, so this control changes the selected-lambda final refit without silently
changing the CV search budget.

## Performance audit

The `p = 100`, C(3), `K = 4` profiles show that essentially all RoCE wall
time is spent in the three nuisance cross-validation paths (initial density
ratio, calibrated density ratio, and calibrated outcome); correction,
aggregation, and variance assembly take milliseconds. The production path
therefore parallelizes independent sources and treatment arms, while keeping
BLAS single-threaded to avoid oversubscription. Federated-DR and Pooled-DR also
reuse one source-to-target density-ratio fit rather than fitting it separately
for each method and arm.

For the binary FACE DGP, the calibration and fixed superpopulation truth use
the same deterministic 100,000-draw target reference population. The generator
computes those moments once per task and reuses them, avoiding a duplicate
`100000 x p` reference draw without changing the RNG stream or estimand.

A 50-point nuisance grid was materially faster in profiling, but it changes
the tuning search and is not a purely computational optimization. The formal
manuscript manifests therefore remain locked to 100 while the five-seed,
provenance-isolated 50-versus-100 study above is incomplete. They may be
regenerated at 50 only if every predeclared numerical, SE, runtime, and
implementation-QC criterion passes; otherwise the 100-point lock remains.

The final 100-point implementation also stops an individual fold's descending
penalty path after three consecutive candidates fail the complete-coordinate
convergence check. The remaining, less-regularized tail is left at its
preallocated infinite validation score. An isolated failure therefore has two
later candidates in which to recover, and a failed coefficient vector is
rolled back before the next warm start. Result files report the skipped-fold
count separately from the legacy invalid-fold and invalid-lambda counts.

This acceleration is guarded by `audit_cv_tail_equivalence.R` and
`run_cv_tail_equivalence_audit.sh`. On the locked `p = 100`, C(3), `K = 4`,
`rho = 2.5`, seed-1 task with 100 nuisance penalties and 5,000 bootstrap
draws, it skipped 17,850 fold fits and reduced elapsed time from 6,559 to
3,516 seconds (1.87-fold). All 100 primary scientific/final-fit fields and all
97 sensitivity scientific/final-fit fields were CSV-exact. The only changed
legacy diagnostics were invalid-candidate counts: every increase was required
to be nonnegative and no larger than its matching explicit skip count. The
audit writes `cv_tail_scientific_equivalence_passed.txt` only after those
conditions hold; it deliberately does not describe the mechanically affected
invalid-count fields as bitwise identical.

## Reused cutoff and inference-radius sensitivity

The primary C(3), `K = 4`, `rho = 0/2.5` tasks write a provenance-checked
sidecar after their main 13-method CSV. One fitted nuisance object supplies
`c = 1, 1.5, 2, 2.5, 3` and
`M_tau_inference = 4, 5, 6, Inf`; correction/IF components are refreshed once
per inference radius, while all cutoffs reuse them. The default `(c,M_inf) =
(2,5)` sidecar row must reproduce the primary TATE estimate and SE within
`1e-12`, and every sidecar records zero nuisance refits.

The 1,000-task cutoff manifest and submitter remain as a standalone fallback:

```bash
scripts/slurm/submit_cutoff_diagnostic.sh
```

Pass a third argument to override the CPU count. With fewer than 8 CPUs the
same estimator is used, but the two treatment arms run sequentially.

Do not submit that fallback after the primary reused sidecars pass their audit;
doing so would duplicate 1,000 expensive nuisance fits.

### Seed-grouped cutoff pilot

For a bounded finite-sample cutoff diagnosis, the seed-grouped workflow fits
one nuisance object per requested rho and reaggregates all requested cutoffs
without nuisance or influence-function refits. It writes to
`grouped_cutoff_pilot/`, never to the primary `raw/` directory, and requires
the tested production package library. First prepare one seed-1 integration
task whose package-primary-cutoff row must reproduce the separately committed production
row:

```bash
sbatch scripts/slurm/prepare_grouped_cutoff_diagnostic.sh \
  results/direct_tate_mc500_b5000/grouped_cutoff_identity \
  C1 2 1 '0;1;2.5' '1;1.5;2;2.5;3' 1
```

Submit and audit that one identity task before preparing the tuning pilot. The
default preparation command creates the C(1), `p = 100`, `K = 2` tuning
manifest with seeds 501--510, which are disjoint from the primary Monte Carlo
seeds 1--500:

```bash
ROCE_BATCH_SIZE=1 ROCE_MAX_CONCURRENT=1 \
  scripts/slurm/submit_grouped_cutoff_diagnostic.sh \
  results/direct_tate_mc500_b5000/grouped_cutoff_identity/manifest.csv \
  results/direct_tate_mc500_b5000/grouped_cutoff_identity/raw
```

```bash
sbatch scripts/slurm/prepare_grouped_cutoff_diagnostic.sh
```

After inspecting the preparation log, submit exactly one tuning seed:

```bash
ROCE_BATCH_SIZE=1 ROCE_MAX_CONCURRENT=1 \
  scripts/slurm/submit_grouped_cutoff_diagnostic.sh
```

Run prefix summaries on a compute node. The identity smoke uses the default
external-production gate, whereas the tuning pilot must name its disjoint
mode explicitly:

```bash
sbatch scripts/slurm/run_grouped_cutoff_diagnostic_summary.sh \
  results/direct_tate_mc500_b5000/grouped_cutoff_identity 1
sbatch scripts/slurm/run_grouped_cutoff_diagnostic_summary.sh \
  results/direct_tate_mc500_b5000/grouped_cutoff_pilot 1 disjoint-pilot
```

Each seed job evaluates `rho = 0, 1, 2.5` and
`c = 1, 1.5, 2, 2.5, 3`. The installed package's primary-cutoff row must be
bitwise identical to the same task's fixed-cutoff fit; the seed-1 integration task additionally checks
the separately committed primary CSV. Summarize the 501--510 tuning pilot with
the explicit `disjoint-pilot` mode. Advance first to five and then ten seeds
only after the preceding prefix passes
`summarize_grouped_cutoff_diagnostic.R`; that mode rejects any seed in 1--500.
The submitter enforces a hard cap of five jobs per call and concurrency two.
Cutoff comparisons are reported with Monte Carlo uncertainty and paired
squared-error differences, rather than selecting a primary cutoff solely from
the observed coverage.

## Truncation sensitivity

Fitting-radius changes require new nuisance fits. Generate the targeted C(3),
`p = 100`, `K = 4`, `rho = 0/2.5` manifest explicitly:

```bash
Rscript scripts/slurm/build_direct_tate_manifest.R \
  results/direct_tate_mc500_b5000/manifest_truncation_diagnostic.csv \
  truncation 500
```

This manifest contains only `(M_fit, M_inf) = (4,5)` and `(6,5)`: 2,000 tasks
across the two rho endpoints. The four `M_fit=5` settings come from the main
task sidecars. Comparison estimators are omitted because they do not depend on
the truncation radius. Submit through the dedicated bounded wrapper:

```bash
scripts/slurm/submit_fitting_radius_diagnostic.sh
```

As with the primary study, advance only after checking each completed batch;
the manifest size is not authorization to submit all rows at once.

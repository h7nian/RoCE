# RHC current-method application, 2026-09-25

The user authorized submission of the real-data analysis. This is a separate
current-method application; all historical analyses remain unchanged.

## Cohort and estimand

Reproduce the previously reviewed insurance cohort: 5039 records, Private target
(1698), Medicare (1458), Private & Medicare (1236), Medicaid (647). The historical
exclusions are No insurance and Medicare & Medicaid. These insurance strata
emulate sites, not different hospitals. K=3 sources and four total strata.
The treatment is RHC and the outcome is 30-day mortality. The contrast is the
target population's treated minus control outcome mean, conditional on the
observational identification and transport assumptions.

The [public data description](https://hbiostat.org/data/repo/rhc) defines RHC on
the first SUPPORT-eligible day; the [catalogue](https://hbiostat.org/data/repo/crhc)
labels clinical measurements from day 1. These descriptions do not establish
that every measurement precedes treatment. Preserve the historical variables
for comparability and record this timing limitation; do not present the analysis
as proof of a causal clinical conclusion.

The shared 61-column OR/weight feature map, pooled covariate-only median/mode
imputation and centering/scaling are the historical specification. They are not
fold-local preprocessing and do not establish an end-to-end private protocol.
The previously observed seven source/arm moment-support flags and substantial
urine-output missingness remain limitations; no covariate or patient is removed
in response to a fitted treatment effect. A clinical timing or preprocessing
revision would require a separately named data specification and matched reruns.

## Frozen primary profile and comparisons

- Two-layer, one-round, source score-derivative calibration, Hou target anchor,
  joint TATE with separate mu1/mu0 weights, cutoff2 (lambda.5).
- Ten outer folds,100 nuisance penalties, proximal Newton, tolerance1e-10 and
  the verified CV certificate. Fitting and inference radius5 preserve the
  historical real-data choice; radius12 is a named sensitivity, not a hidden
  substitution of the simulation setting.
- Ordinary target-only uses the same target outer folds as RoCE. The calibrated
  anchor has a distinct label. SS, IVW, Federated-DR and Pooled-DR use the same
  records and preprocessing with their existing estimation procedures and
  5000 multiplier draws. Baselines retain their own nuisance-fold procedures;
  identical nuisance folds across all methods are not asserted.
- Each baseline has its own job and saved stages. Check both-arm contrast and
  covariance identities, and later match the ordinary-target reference and data
  hash to the primary RoCE run. Multipliers are not full-pipeline bootstrap fits.

## Predefined supplementary profiles

Eight method profiles: primary; nuisance CV1se; radius12; standard source with
calibrated target; ordinary target with calibrated sources; outer fold seeds
101,202,303. The primary preserves the historical sample-size-based outer folds.
The existing `seed` only initializes the loader; `fold_seed` explicitly changes
outer partitions. All profiles keep the same records. Repeated folds are
stability checks on one cohort, not independent repetitions with known truth.

Primary nuisance fits also support reaggregation at cutoffs1/2/3 under joint,
separate-arm and common weights. No nuisance refit is needed for this diagnostic.
Do not select a primary setting by effect size, statistical significance or
the smallest reported standard error. Failures are retained; there is no
automatic min-to-1se numerical fallback.

## Execution and outputs

Twelve individual Slurm jobs: eight method profiles and four baselines. Four
CPUs,24GB and6h per job; preempt/public/saffo-2tb routes, requeue and exact
nuisance/stage checkpoints. The checked simulation estimation library remains
immutable; the frozen RHC interface is loaded in a separate environment whose
parent is that namespace. Interface tests and every input/workflow/library hash
are recorded. The existing job ownership and duplicate-submission guard is reused.

Each task retains methods, arm-specific source weights, two-arm covariance,
fit/clip diagnostics, warning/error records and checkpointed fitted objects.
Review estimates, CIs, borrowing and numerical/support behavior; do not report
RMSE or coverage without a known population effect. Completed results are not
automatically substituted into the manuscript.

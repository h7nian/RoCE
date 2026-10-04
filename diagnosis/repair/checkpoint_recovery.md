# Preemptible repeat jobs

The low-dimensional v5 recovery described below completed all 1800 repeats.
The user's current high-dimensional batch is
`implementation/r7/validation_highdim_mc200_cutoff2_v1`, with cutoff 2
(`aggregation_lambda=0.5`). The prior cutoff-1 workers were stopped before
their immutable completed-fit checkpoints were staged in `checkpoint_seeds/`.
Each new worker restores only these RDS entries into its own checkpoint
directory; existing restart progress is retained. Input, installed-code and
model hashes are checked by the R cache before a fit is reused. No old
treatment-effect estimate is imported across cutoff settings. Source entries
remain available after the new worker removes its own checkpoint links.
Fresh versus restored cutoff-2 results were checked to 1e-12, including both
arms running concurrently and all three aggregation modes.

Slurm copies batch scripts into its spool directory, so the controller cannot
find neighboring Python code through `dirname "$0"`. The wrapper now accepts
the frozen workflow directory explicitly, defaulting to the campaign workflow.
An integration test runs the script from a separate spool-like path, including
paths with spaces. Worker estimator files from earlier campaigns are unchanged.

## Historical low-dimensional recovery

The active recovery campaign is
`/scratch.global/zhan9381/FACE-HD/implementation/r7/validation_mc200_v5/`.
It references 1628 completed v4 repeats with checksummed import records;
only the remaining 172 original seeds require computation. The v4 logs show
node startup failures and filesystem I/O errors on acl03, acl41, acl74,
acl83 and acl96. The last node entered DRAIN after a health-check watchdog
failure. `--exclude-nodes` excludes these explicit node names from both
workers and the controller while retaining the requested partition list.
Several old failed job IDs are no longer available to `scontrol requeue`.
The recovery therefore uses new task directories/job IDs for those seeds;
their original partial checkpoints and errors remain in v4. This is a fresh
same-seed computation, distinct from the normal same-job checkpoint resume
described below. Scientific settings, installed library and worker R/shell
contents are unchanged. The original controller 1497946 was cancelled.

The campaign contains 1800 planned repeats: p10/20/50 x C1/C2/C3 x 200,
K2, 1000 observations/site. Imported completed results retain their original
job records and file hashes and are not recomputed.
The user subsequently re-enabled `saffo-2tb` and requested this routing list:
`preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb`.
Partition changes within a campaign use the optional `routing.json`
operational override and update pending jobs without changing estimator
inputs or checkpoint fingerprints. The listed order is not a scheduling
priority. New preparation defaults use the same list.
During the earlier v4 rollout, the same monitor ran on the
already allocated ahn1 node, with Slurm job 1497946 as a fallback. Its bounded
local lease and allocation are recorded in `local_controller.json`. If the
allocation ends, the fallback dependency permits takeover; on an ordinary local
lease expiry the launcher releases that dependency. No new estimator jobs run
inside the local monitor.

## What a checkpoint contains

`checkpoint_dir` reuses the existing exact-input nuisance cache in
`run_single_simulation()`, `run_tate_crossfit()` and `run_crossfit()`.
Each completed nuisance fit is atomically saved with its input/model checksums.
The persistent directory is also keyed by the installed R/native-code hashes.
A restarted repeat regenerates the same data and folds from the same seed,
recovers matching completed fits, and recomputes any unfinished fit and the
remaining assembly. It does not save an optimizer in the middle of an
iteration or claim to resume an entire R process image.

The argument requires `use_lambda_cache=TRUE` and cannot be combined with
`nuisance_cache_dir`, whose older temporary-cache behavior remains unchanged.
Default `checkpoint_dir=NULL` preserves the ordinary fitting interface.
Completed result files are committed before intermediate checkpoints are
removed. First-seed full fitted artifacts remain available for diagnostics.

## Slurm and submission state

Workers use `--requeue`, `--signal=B:USR1@120`, appended logs, and the same
Slurm job ID on restart. Task claims check job identity, configuration hash,
and restart count; a worker lock prevents concurrent writers. Completed
repeats exit without rerunning. Genuine model errors are not silently retried.
The default maximum is 20 restarts. Prior attempt errors/status files are
retained under each task's `attempts/` directory.

The controller records each job submission and recovers original job IDs by
unique scheduler comments when a submission was interrupted. Ambiguous
submissions require review rather than duplicate submission. The live queue
supersedes stale PREEMPTED accounting states. A user cancellation stops the
controller; it does not deliberately requeue itself after SIGTERM.
The configured one-hour controller requests requeue at its safe time budget,
or after its advance USR1 signal. Batch size and scheduler request pacing are
separate from the 512-job in-flight limit.

MSI's current `preempt` policy has `PreemptMode=REQUEUE` and `GraceTime=0`.
Checkpoints are therefore written during normal fitting, rather than relying
on a warning before preemption. Slurm's advance `--signal` is associated with
the requested time limit; it does not by itself guarantee an early preemption
warning. See [Slurm sbatch](https://slurm.schedmd.com/sbatch.html).

The first controller encountered a metadata-visibility race: Slurm reported
COMPLETED before a shared-filesystem completion marker became visible.
Controller workflow v2 now waits up to 120 seconds before classifying a missing
marker as an error. Tests cover both delayed visibility and a genuinely
missing result. Controller 1489557 replaced 1486262; the worker jobs and their
frozen implementation were preserved.

## Validation

The new checkpoint tests and related cache/three-level tests passed 10 cases
and 140 assertions, including a deliberate interruption and restoration to
within 1e-12 for estimates, standard errors and the three aggregation modes.
Sequential and concurrent-arm restoration were checked. Workflow tests cover
job ownership, duplicate prevention, accounting/requeue transitions and
recovery of a previously submitted job ID.

The first two direct-preempt repeats completed normally. The p10/C3 repeat
matches the previous non-checkpoint result exactly in the saved estimate and
standard error. These first jobs had no actual Slurm preemption; the numerical
interruption test is recorded separately from scheduler validation.

The earlier v2 jobs 1482279 and 1482280 failed on acl03 during startup;
1482279 logged `TaskProlog failed status=1`. Their records are retained.
Neither had produced an R estimate; the replacement repeats completed under
1485127 and 1485128. No statistical failure is discarded from a performance
summary, and one-repeat diagnostics do not establish Monte Carlo coverage.

## Queue interpretation

`agsmall` and `amdsmall` remain enabled legacy partitions but share hardware
with the broader CPU partitions. Pending counts and one representative node
are not a sufficient estimate of scheduling delay. A comma-separated
partition list is not a priority order, and a short list does not establish
that Slurm was previously pinned to one partition. Sources:
[MSI shared partitions](https://userdocs.msi.umn.edu/compute/shared_partitions.html),
[MSI March 2025 bulletin](https://msi.umn.edu/getting-started/help/users-bulletin/users-bulletin/msi-users-bulletin-march-2025),
[Slurm sbatch](https://slurm.schedmd.com/sbatch.html).

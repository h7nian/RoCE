# Current-paper fixed-K MC200 expansion

The user selected two-layer cross-fitting for this expansion after the completed
paired layer comparison found no systematic TATE advantage from three layers.
This is an operational choice supported by those experiments, not an equivalence
theorem. The existing three-layer implementation and historical results remain
available for comparison and next-paper growing-K research.

## Primary panel

| Control | Prespecified setting |
|---|---|
| Observations | 1000 at every site |
| Source count K | 2,4,6,8, in addition to one target |
| Covariate dimension p | 100,200 |
| Scenarios | C1,C2,C3 |
| Source outcome shift rho | 0,.5,1,1.5,2 |
| Shifted sources | Source1's treated arm, when rho>0 |
| Repeats | Seeds1–200 in every cell |
| DGP | Frozen bounded_joint_v3, four active slopes, fixed population truth |
| Training | Two layers, ten folds, fold-summed calibration, one-round protocol |
| Target anchor | Hou calibrated |
| Primary aggregation | joint_tate; separate_arms/common_tate retained |
| Aggregation cutoff | 2, with lambda=.5 |
| Nuisance fitting | 100 lambdas, proximal Newton final tolerance1e-10, radius12 |

There are120 primary cells and24000 two-layer repeats. Each is paired with
ordinary target-only, sample-size weighting, inverse-variance weighting,
federated DR and pooled DR on the same generated data. The calibrated target
anchor is also retained. Complete old baseline results are reused only for
matching scientific/data identities; final reviews must still verify data
hashes and target estimates/SEs/truths. Old three-layer estimates cannot count
as two-layer repeats, and the earlier401–450 layer seeds are not relabeled1–200.

The checked CV certificate is explicitly enabled in new runs. It retains the
full penalty grid and numerical acceptance rule; full-path, pipeline and baseline
parity checks are archived under r10. Coordinate descent remains available.

## Supplementary panels

- C4, both nuisances misspecified: p100/200, K2/4/6/8, all five rho values,
  200 repeats/cell, with matching baselines. This is outside the nuisance
  correctness union and is labeled as a boundary.
- Half-invalid sources: p100/200, C1–C3, K4/6/8, rho.5/1/1.5/2,
  200 repeats/cell, with matching baselines. K2/one-invalid and rho0 controls
  already belong to the main panel and are not resubmitted.
- Source final-calibration ablation: K2, p100/200, C1–C3, rho0,
  200 paired repeats, with the Hou target anchor held fixed. The standard source
  retains initial merged-weight balancing and ordinary source OR fitting.
- Target-anchor ablation: the same six cells and200 repeats, with calibrated
  sources and ordinary target fitting. Both ablations retain outer-training eta
  learning, cutoff2, nuisance settings and generated data.

The two single-component ablations are meaningful two-layer comparisons. The
historical four-cell calibration factorial remains separately labeled: removing
final calibration from both source and target eliminates the calibration role
itself, and is not silently treated as the same two-layer program.

The combined plans contain46400 method/baseline data settings, plus1200 source
and1200 target ablation results. Completion inventory and reuse counts are in
the frozen plan; these counts are not independent Monte Carlo draws across
methods, rho values or scenarios. Different K changes total sample size and
source tilt positions. Cross-K seed labels do not imply identical datasets.

## Execution and review

`prepare_current_paper_mc200.py` reuses the checked launcher and coordinator.
One repeat remains one Slurm job. New routing includes
`preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb`, with
checkpoint/requeue and the documented failed-node exclusions. Partition order
does not impose scheduler priority. Main panels have higher concurrency;
supplements wait for the matching primary pipeline's first-seed execution
checks and use smaller concurrency caps. Checks never accept/reject seeds on
coverage, bias or RMSE. Existing submitted work is not duplicated.

The controller caps are summed and recorded, below the verified5000-job account
limit with room for other work. Scientific configurations, installed libraries
and worker scripts are frozen before submission. All results, logs and snapshots
are under `/scratch.global/zhan9381/FACE-HD/implementation/r11/`.

Review all prescribed seeds, paired data, target references, actual source
counts, mean/variance identities, nuisance convergence and the three aggregation
modes. Report bias, RMSE, empirical SD, mean SE, coverage with Monte Carlo
uncertainty, interval length and actual arm-specific source weights. Weak
departure is a boundary of the current method, not a reason to tune the DGP,
filter repetitions or require every mode to perform well.

## Calibration mechanism and next-paper scope

Use the saved outer/inner score records to compare source candidates before
aggregation, both arm means and their covariance, and calibrated target-only.
Add comparisons with the same outer-training-derived eta applied to each
candidate program, recomputing each method's variance. Evaluate score derivatives
on independent diagnostic observations for the saved fitted models; small
gradients on the very calibration training rows would be a circular check.
These analyses distinguish nuisance-error control, variance changes and the
effect of relearning aggregation weights.

Fixed-cutoff sensitivity can reuse sufficiently complete saved fits or score
records, with eta and its SE terms recomputed. A sample-size-dependent cutoff
study is separate from this1000/site panel. No cutoff is selected from the
same reported Monte Carlo outcomes.

The next paper must address weak-departure inference and growing K with explicit
uniform nuisance/covariance/bias conditions and an efficiency analysis. Three
layers are available for that work but are not themselves a solution. The
current paper need not solve those open problems before its experiments are
reported honestly.

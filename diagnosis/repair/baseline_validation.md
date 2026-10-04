# Baseline validation and paired runs

The existing cutoff-2 RoCE run contains ordinary target-only AIPW and the
calibrated target anchor. The supplementary baseline run adds the four
comparisons already described in the manuscript: sample-size aggregation,
inverse-variance aggregation, density-adjusted federated AIPW and pooled
density-weighted AIPW. It uses the same frozen RoCE library, generated records,
site sizes, working features and seeds. The target reference uses the identical
target folds and ordinary nuisance fits used by the main run.

The nine pipeline checks cover C1/C2/C3 at p10, p100 and p200, with 1000/site.
All completed. Each checks finite estimates/SEs, the treatment-minus-control
identity, direct reconstruction from site influence blocks, and the two-arm
variance/covariance identity. Available main-study data hashes, target estimates
and target SEs must match to 1e-12. Package regression checks independently refit
case-weighted nuisance models, including density normalization/clipping, and
check bootstrap reproducibility and covariance handling. Warnings from small
test fixtures are retained in the check logs.

The baseline worker checkpoints the target reference, density fits, arm-1
results and TATE results separately. Each checkpoint includes a data/config/task
hash, a value checksum, warnings and the post-stage RNG state. A restart reuses
completed stages and restores the RNG before unfinished work. It does not
rerun the expensive RoCE calibration merely to obtain baselines.

## Interpretation of the implemented comparisons

| Comparison | Point estimator | Variance used for TATE |
| --- | --- | --- |
| Target-only | Held-out ordinary AIPW on target records | Paired target influence scores |
| SS | Sample-size average of local AIPW arm means, then difference | Site-aligned paired influence blocks |
| IVW | Local AIPW arm means with fitted DL-adjusted inverse-variance weights | Paired influence blocks conditional on fitted meta-analysis weights |
| Federated-DR | Density-weighted local AIPW, combined with fitted IVW weights | Shared-target and source influence blocks with nuisance/density corrections |
| Pooled-DR | Pooled nuisance fits and density-weighted AIPW | Site-stratified paired influence blocks with nuisance/density corrections |

Comparison intervals use 5000 multiplier draws, retaining their analytic
influence variances as diagnostics. This is an influence-function multiplier
bootstrap, not a bootstrap that refits every nuisance and selection decision.
The TATE assembly uses both arms jointly; it does not add marginal variances
or add between-site heterogeneity to conceal target-population bias.

SS/IVW local means generally describe different populations under covariate
shift. The implemented federated/pooled weighting benchmarks rely on a suitable
density adjustment; their labels do not assert equivalence to every published
federated or fully transport-doubly-robust estimator. In the bounded joint DGP,
a correct merged arm/site tilt does not imply that the separate marginal
source-to-target density ratio is log-linear. Model-class limitations and
conditional active-set/weight inference remain distinct from coding errors.
All performance, warnings and failed repeats are retained for evaluation.

The paired full run is
`/scratch.global/zhan9381/FACE-HD/implementation/r7/baseline_mc200_v1/`:
p10/20/50/100/200 x C1/C2/C3 x 200 seeds. The original low-dimensional main run
used cutoff 1 and the high-dimensional main run uses cutoff 2; these cutoffs
do not enter the four baseline estimators. Pilot and regression outputs are
in `baseline_pilot_v2/` and `baseline_checks_v1/` under the same r7 directory.

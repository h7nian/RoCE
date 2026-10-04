# Two/three-layer comparison

The two-layer definition follows the outer-evaluation / secondary-calibration
structure in the reviewed Overleaf `JASA/JASA.tex` (revision `d57e90e5`). The
specific eta step is completed as described in `two_layer_jasa_review.md`:
reuse final outer nuisance fits on the outer training sample. These are
training scores, not independent validation scores. Each outer fold still
has its own eta and is excluded from every statistic used to learn it.

`crossfit_layers=2L` selects that program. `crossfit_layers=3L` independently
validates complete calibrated inner fits. Both use the same outer estimator;
the old `initial` source-validation ablation is a different program. Existing
production campaigns and package compatibility defaults remain unchanged.

## Initial paired design

| Control | Prespecified values |
| --- | --- |
| Sample size | 1000 at each site |
| Original folds | 10 |
| Covariate dimension p | 10, 100, 200 |
| Source count K | 2, 4, 6, varied separately from p |
| Scenarios | C1, C2, C3 of bounded_joint_v3 |
| Outcome deviation rho | 0 in this initial layer comparison |
| Repeats | 50 per cell and layer, seeds 401–450 |
| Nuisance program | Source and target score-derivative calibration |
| Communication | One round |
| Nuisance optimizer | Proximal Newton, tolerance 1e-10, grid size 100 |
| Aggregation cutoff | 2 (lambda 0.5) |
| Aggregation modes | joint_tate primary; separate_arms and common_tate comparisons |

There are 27 data cells and two layer choices: 2700 prescribed repeat results.
Each worker evaluates all aggregation modes using the same nuisance fits.
The first local pair is imported into the p10 campaigns, not rerun. Matched
ordinary and calibrated target-only estimates come from each worker. The
existing current-paper baseline, deviation and source-count campaigns continue
separately; this layer comparison does not relabel their results.

Execution proceeds from a complete ten-fold, 100-lambda local pair to the p10
first-seed wave, then to the remaining p10 and high-dimensional repeats. Gates
check complete numerical outputs, matching data and final outer-fit fingerprints,
and complete arm variance diagnostics. Coverage and RMSE never determine which
seeds are released, kept or discarded. Each repeat is one Slurm job with the
existing checkpoint/requeue protocol and permitted partitions.

## Estimation and inference checks

Retain the actual aggregated mu1, mu0 and TATE for all three aggregation modes,
plus both arm variances and their covariance. The arm contrast must reproduce
the reported TATE estimate and variance, both with fixed weights and with the
implemented eta sensitivity. Tests compare the latter with empirical-mass
finite differences that relearn eta while keeping fitted nuisances fixed.

Reports include bias, RMSE, empirical variance across repeats, mean estimated
variance, variance ratios, coverage intervals, paired MSE differences with
Monte Carlo standard errors, and eta distributions. The distributions of
repeat-mean weights are kept separate from the fraction of outer-fold weights
numerically near zero. Runtime comparisons are descriptive unless hardware,
parallelism, restart and cache conditions match.

Fifty repeats are an initial diagnostic, not enough to establish small coverage
or variance differences precisely. Further repetitions should improve precision
of the same prespecified comparisons. Deviations, weaker separation and growing
K need additional experiments and theory; this initial rho-zero design does not
establish their validity.

## Completed comparison

All 2,700 prescribed results are complete: 1,350 matched two-/three-layer
pairs, with identical data and final outer-fit fingerprints. The final combined
tables are under scratch `implementation/r9/layer_review_p200_complete_v1/`.
Its `p200_interpretation.md` and `figures_p200/` describe the last completed
panel; the tables also retain p10 and p100.

For joint TATE, every paired MSE contrast in all 27 cells is within two Monte
Carlo SEs. This is evidence against a large, systematic advantage of the third
layer in these fixed-K, rho-zero settings, not an equivalence test. At p200,
the maximum relative RMSE difference is1.25%, and coverage is identical in
seven of nine cells and differs by one repeat in the other two.

Mean-specific and variance results are not all negligible. At p200/C3/K4,
two-layer mu0 RMSE is2.64% lower (paired MSE contrast3.70 MCSEs), while
three-layer reported variance is closer to empirical variance. For TATE,
p200/C1/K2 and C3/K4 have paired variance-discrepancy contrasts about2.47 and
2.65 jackknife SEs; both layers overestimate in those cells and three is less
conservative. These exploratory contrasts are not multiplicity adjusted.
Three layers do not systematically correct the observed undercoverage.

The experiment supports retaining two layers as a serious fixed-K option.
The training-moment consistency argument below is still needed, and the
existing production default remains unchanged. Recorded runtimes include
different machines, checkpoint reuse and restarts, so their ratios do not
measure the isolated computational cost of an extra layer.

## Theory scope

For fixed K, consistent eta and valid-source root-n differences can make eta
estimation a higher-order contribution after incompatible sources are removed.
Two-layer reuse still needs consistency of its actual training-score means and
covariances. Independence of the final outer fold does not prove that lemma.
Fixed cutoff2 is also distinct from the manuscript's growing-cutoff oracle
selection theorem. The reported eta sensitivity is conditional on nuisance fits
and active sets; empirical repeated-simulation variance is the primary check.

The current lower-level package reports the realized layer count as
`crossfit_levels`; public arguments and the new diagnostic CSVs use
`crossfit_layers`. This preserves existing result readers without introducing
another competing algorithm option.

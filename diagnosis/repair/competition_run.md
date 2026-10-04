# High-dimensional competition results

The user selected a fixed aggregation cutoff of 2 on 2026-09-20, after the
low-dimensional cutoff-1 review. Later authorization includes baselines and
conditional deviation/ablation validation after execution checks pass.
The active specification is p100/p200 x C1/C2/C3 x 200 original seeds, with
1000 observations at the target and at each of two source sites. The bounded
joint DGP, sparsity 4, covariate/treatment shift controls of 1, rho 0, ten folds,
three-level fitting and one-round protocol remain unchanged.

`--aggregation-cutoff 2` sets `aggregation_lambda=0.5`. For discrepancy Wald
statistic T, the source coefficient penalty is `(0.5 * T - 1)_+ * abs(eta)`.
This activates soft shrinkage beyond two estimated standard errors; it is
neither a finite-degrees-of-freedom t test nor an automatic source exclusion.
The 100-point nuisance lambda grid and solver tolerance 1e-10 are separate
settings and remain fixed.

Prioritize joint-TATE results, compare with separate-arm aggregation and the
target-only benchmark, and retain common-TATE as the constrained comparison.
Report every planned repeat, any execution failures, coverage with Monte Carlo
uncertainty, bias, RMSE, empirical SD, mean SE and their ratio. Compare errors
on matched seeds. Do not choose or discard repeats based on their estimates.

The new results belong to
`/scratch.global/zhan9381/FACE-HD/implementation/r7/validation_highdim_mc200_cutoff2_v1/`.
The stopped cutoff-1 campaign and its outputs remain separate. Only immutable
completed nuisance fits can be reused across the cutoff change. Their exact
inputs, installed-code namespace and model checksums remain checked by R.
Workflow checks passed 13 Python cases; two R checkpoint cases passed 14
assertions, including cutoff-2 cached versus fresh estimates and standard
errors within 1e-12 for all three modes and concurrent arms.

The completed low-dimensional results used cutoff 1. Label this difference
explicitly if those results are shown alongside the new study. Current rho-0
experiments do not evaluate nontransportable outcome mechanisms. Fixed-cutoff
results also do not establish the growing-cutoff oracle theorem in the saved
manuscript; joint-TATE theory needs its own corresponding conditions. Additional
settings are kept in separate paired campaigns under
`implementation/r7/extended_validation_v1/`, retaining cutoff 2. Their gates
check execution and data integrity; favorable coverage is not a release condition.

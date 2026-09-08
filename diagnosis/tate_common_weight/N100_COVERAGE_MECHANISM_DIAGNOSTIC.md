# Saved n=100 coverage-mechanism diagnostic

This is a post-checkpoint diagnostic, not a new simulation or a CI-selection
gate. Inputs are the frozen 100 v19 seed bundles and their n=100 summary.
No seed, observation, estimator, standard error or published interval is changed.

For each saved primary TATE fit, independently reconstruct target-only plus
each source's exact fold-specific contribution using the original observation
indices and site sample sizes. Report source weights separately from these
contributions; never substitute mean weights times pooled source estimates.
Verify the sum against both the saved fit and raw CSV. Across the same 100
paired seeds, report bias and variance decompositions of target plus borrowing.
These algebraic decompositions do not by themselves prove why a source is biased.

Separately report actual interval misses above/below truth, studentized errors,
and diagnostic fractions after subtracting the observed Monte Carlo mean bias
or replacing standard errors by the observed Monte Carlo SD. These use truth
and the same 100 outcomes: they are deliberately oracle/descriptive thought
experiments, NOT valid confidence intervals, bias corrections, or coverage
validation. Do not choose a future method, cutoff or inflation factor to make
these fractions look favorable. They only separate observed location, scale
and tail-shape symptoms before a separately justified inference revision.

The exact variance decomposition is Var(target + borrowing) = Var(target) +
Var(borrowing) + 2 Cov(target, borrowing). It is not an independent additive
weight-variance correction. Bias decomposes additively into target error and
the two source increments. Sampling MCSEs use independent seed-level values,
never folds or paired rho settings as additional independent replications.

## Completed saved-data diagnostic

Job 18410137 completed 0:0 in 2:26. The checksummed output is
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/coverage_mechanism_n100_v1/`.
All 600 target-plus-source reconstructions match saved fits/CSVs, maximum error
4.996e-16. All 6,000 fold/source contributions are retained. Bias/variance
decompositions close to numerical precision. No nuisance fit was rerun.

| Rho | Target mean error | Mean s1 increment | Mean s2 increment | RoCE mean error | Mean s1 weight |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | -0.011766 | +0.000100 | -0.000351 | -0.012017 | 0.31387 |
| 0.5 | -0.011766 | +0.022851 | -0.000303 | +0.010782 | 0.26417 |
| 1 | -0.011766 | +0.022767 | -0.000427 | +0.010575 | 0.14934 |
| 1.5 | -0.011766 | +0.011929 | -0.000586 | -0.000423 | 0.05988 |
| 2 | -0.011766 | +0.005792 | -0.000635 | -0.006608 | 0.02508 |
| 2.5 | -0.011766 | +0.002496 | -0.000624 | -0.009894 | 0.00935 |

At moderate rho, the perturbed s1 still receives nonzero weight and shifts
the estimator upward enough to reverse the negative target-only error.
This identifies the numerical source of the observed bias, not a proof that
every increment estimates a population bias or that a different cutoff is
justified. Weights are only descriptive averages; the increments use actual
fold/site fractions and are not mean weights times pooled estimates.

At rho=1, original intervals miss below/above truth through estimates that
are too low/high in 3/15 seeds; the relearned diagnostic misses in 1/11 seeds.
The mean analytic variance is 0.00055244 versus repeated-sampling variance
0.00092865; mean relearn variance is 0.00088109. At rho=1.5, these variances
are 0.00063577, 0.00111024 and 0.00099623, respectively. Mean analytic
variance is therefore only about 59%/57% of empirical variance in these settings.

The truth-assisted thought experiments further separate symptoms:

| Rho | Original fraction covered | Remove MC mean bias, keep analytic SE | Use MC SD, retain bias | Remove bias and use MC SD |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.82 | 0.85 | 0.94 | 0.96 |
| 1.5 | 0.83 | 0.83 | 0.95 | 0.95 |

These are NOT corrected coverage estimates or deployable intervals: the
centers/scales use the same 100 errors and known truth. They show that bias
removal alone cannot explain the original shortfall, especially at rho=1.5.
Removing the MC mean bias while retaining relearned SEs gives diagnostic
fractions 0.93/0.92. The studentized relearn errors have SD 1.141/1.126,
despite the closer mean SE scale. Matching average variance is insufficient
to validate the distribution of each estimate divided by its own SE.

The next inference analysis must handle both retained source discrepancy and
the joint uncertainty of data-adaptive weights and evaluation scores. The
identity Var(target+borrowing) includes a negative covariance term here;
simply adding an independent variance correction is not supported.

General theory explains why successful computation alone is not enough:
orthogonal-score/cross-fitting inference still requires nuisance-rate and
sampling conditions ([Chernozhukov et al.](https://doi.org/10.1111/ectj.12097)),
and an ordinary bootstrap is not automatically valid for nonsmooth maps
([Fang and Santos](https://arxiv.org/abs/1404.3763)). These sources are general
requirements, not a theorem already established for this RoCE implementation.

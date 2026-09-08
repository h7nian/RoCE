# Independent fixed-fit population check

Hold the selected C2/C3 outer-1, source-s1, rho-0 nuisance and projection
coefficients fixed. Generate 25 independent batches of 2,000 new covariates
per site with the frozen FACE K=2, p=100 generator, using seeds 918731+b
for batch b=1,...,25. Generate both source sites to retain the K=2 skewness
definition, but evaluate target and source s1 only.

Integrate A and Y analytically conditional on X: replace each arm indicator
by its true conditional probability and Y in arm-weighted residuals by that
arm's true conditional mean. The score and derivative-row expressions are
affine in A and A*Y, so this computes their conditional expectations without
Bernoulli draws. Verify that identity before running the integration.

Retain batch-level score means, target conditional truth, and coordinate
gradient means. This estimates the population behavior of these fixed fits,
not repeated-training bias or coverage. It uses 50,000 new X values per
evaluated site, no bootstrap and no nuisance refits. Do not change selected
scales or clipping after seeing its results. Integration uncertainty is not
a standard error for the fitted estimator.

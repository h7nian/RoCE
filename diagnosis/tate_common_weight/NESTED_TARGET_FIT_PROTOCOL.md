# Three-fold-excluded target nuisance fits: fixed development policy

Use the frozen v19 seed 10013 target data and outer treatment-balanced folds.
For each sorted excluded triple (a,b,c), train on the remaining two folds.
Share exactly one P(A=1) fit and two arm-specific binary outcome fits across
all requests with that excluded triple. Do not reuse a one-/two-fold-excluded
model. This creates 30 fits for 60 ordered validation requests.

The model seed is 730000 + 100*(100*a + 10*b + c) + model_index, with indices
1=propensity, 2=treated outcome, 3=control outcome. Generate explicit balanced
CV fold labels by shuffling repeated 1:n_cv, using the package CV fold count
with min_per_fold=10. Fit binomial lasso, 100 lambdas, package maximum
iterations, lambda.min, and no parallel CV. This is a new deterministic fit
policy, not a replay of nonexistent saved three-fold-excluded models.

Retain full CV fits, coefficient vectors, observation IDs, CV labels, seeds,
all warning messages and failures. Insufficient binary class support in any
CV training split fails that model. Nonzero glmnet jerr or failed prediction
reconstruction prevents downstream eligibility. Warnings remain visible even
when a fit passes these numerical checks. Do not change seeds, choose another
lambda rule or omit failed subsets after inspecting results.

This bounded fit gate does not select the final estimator, establish nuisance
rates, or validate inference. Inner projection candidate selection will use
these fits only after the exclusion and numerical checks pass.

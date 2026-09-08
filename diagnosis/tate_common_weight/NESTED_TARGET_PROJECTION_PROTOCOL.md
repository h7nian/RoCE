# Inner target projection validation policy

Apply the existing target projection risk, candidate scales {0.5,1,2}, null
candidate, coordinate SD penalties, numerical tie rule and all-fold failure
policy to each of 20 ordered outer/evaluation fold pairs. Validate over the
remaining three folds using the audited three-fold-excluded nuisance fits.
Fit projection coefficients only on those models' two-fold training subsets.
Retain all 240 candidate/validation attempts, including failures.

After selection, fit final projection coefficients on the three-fold training
complement using the original cached two-fold-excluded nuisance models.
Only then inspect the inner evaluation fold. Check its predictions and
unadjusted TATE against the original saved inner fits. Retain training and
evaluation clipping counts, selected fits, every observation score and IDs.

Do not reuse outer-selected tuning parameters, tune against TATE truth,
discard failed pairs, replace clipping derivatives or recompute production
weights in this check. This is exploratory method development on seed 10013;
the source graph and full fitted-estimator inference remain unresolved.

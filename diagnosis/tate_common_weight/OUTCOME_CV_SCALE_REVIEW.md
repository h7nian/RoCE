# Outcome CV/refit normalization discrepancy

Inspection of the calibrated and refined outcome CV kernels shows training
weights multiplied by n_arm_train/n_source_full. Their GLM update divides
by n_arm_train, so the training data term is a partial sum divided by the
original full source size. The final fit uses n_arm/n_source_full times the
arm-average loss. On balanced five-fold CV, the training data term is thus
0.8 times the corresponding source-population-normalized training mean.

At a fixed reported lambda, the CV training objective is equivalent to a
source-population-normalized objective with lambda/0.8, whereas the full fit
uses lambda. Validation losses also use full-source denominators; their
common scale in equal folds does not change min/1SE rankings by itself.
Do not conflate validation-score scaling with training penalty scaling.

`check_outcome_cv_scale.R` verifies the frozen native calibrated CV kernel
against closed-form Gaussian lasso solutions on a centered paired design.
All three native CV scores reconstruct within 1.74e-18, all 15 fold/penalty
equivalences pass, and the final native fit agrees with its full-source
closed form within 1e-7. The centered design isolates scaling from the native
CV tolerance floor (1e-4); an earlier noncentered fixture had score error
2.42e-7 and was not treated as a scaling contradiction or used to relax the
test threshold.

This is a verified CV/refit objective-normalization discrepancy relative to
using the same source-population loss scale. It is **not** an established
explanation for the large selected C3 penalty or coverage loss: the current
full refit is effectively less penalized than its CV training candidate.
Correcting this alone need not improve the observed bias. Actual shared-fit
objects retain selected min/1SE penalties and failure counts, not the complete
CV score paths; their particular rankings have not yet been independently
reconstructed.

The bundle `outcome_cv_scale_audit_v1/` under the independent-pilot root has
three verified payload checksums. No production or frozen v19 source was
modified. Any proposed normalization correction must be tested in an isolated
package candidate across calibrated/refined and binary/continuous paths before
integration, with fresh full R tests and statistical gates.

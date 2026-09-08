# Observation-level derivatives of selected corrected source scores

`source_score_derivative_rows.R` computes the derivative with respect to
every initial/final nuisance coefficient, holding the projection coefficients
fixed. It preserves raw versus clipped tilting terms, the target/source
moment allocation, and original calibration-piece weights. The target and
source row means sum to the complete D - J' a derivative.

On the selected C2/C3 outer-fold scores, the four arm/configuration systems
passed 48 site/direction finite-difference comparisons. Maximum row error
is 5.71e-9 and maximum aggregate derivative identity error is 6.77e-17.
The executable check is `check_source_score_derivative_rows.R`.

For context only, the check also scales each residual gradient by
sqrt(sample_var_target/n_target + sample_var_source/n_source). Maximum
absolute ratios are 3.2596/3.2202 for C2 control/treated and 2.9520/3.6053
for C3 control/treated. These maxima range over 1,510 coordinates. They are
not single-coordinate z tests, simultaneous bounds, or valid full-estimator
SEs. Fitted nuisance dependence and treatment-balanced fold sampling are not
accounted for by this descriptive scale.

Therefore the earlier increase in outer-fold maximum raw gradients is
insufficient to establish population orthogonality failure on its own.
Conversely, noisy gradients cannot justify claiming successful orthogonality.
Population/remainder analysis and independent sampling gates are still needed.
No candidate, penalty, observation or reported interval was changed in response
to these diagnostic ratios. The derivative rows are an ingredient for further
analytic covariance work, not a complete influence function.

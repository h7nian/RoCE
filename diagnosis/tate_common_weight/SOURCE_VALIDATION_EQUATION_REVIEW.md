# Source equation evaluation on validation observations

`source_projection_validation_state.R` constructs an evaluation-only view
of a fitted source equation. It rejects observation overlap with any initial
training or final calibration subset. It does not refit parameters or alter
the coefficient graph.

Initial nuisance moments each use the validation-site empirical mean. Final
target calibration summaries retain equal weights across the original
calibration pieces. Final source calibration moments retain each original
piece's sample fraction. Consequently the same validation rows can evaluate
each plug-in without multiplying the full source moment by the number of
pieces. The source denominator for piece k is n_validation / fraction_k.

The existing equation evaluator accepts these explicit denominators only
when supplied; its default training calculation is unchanged. The regression
check compared both first-split arm systems to their saved evaluations using
`identical`, including Jacobians and gradients, and passed.

`check_source_validation_state.R` then compared validation moment products
against the separately implemented per-observation holdout correction for
six deterministic coefficient directions per arm. Maximum discrepancy is
4.34e-17. Validation directional Jacobian errors are at most 8.34e-10 and
score-gradient errors at most 9.78e-11. Deliberate training-observation reuse
was rejected. This uses the audited first source validation split only.

These checks supply the held-out equation needed for source coefficient-risk
evaluation; they do not select a penalty or prove asymptotic validity. Other
validation fits and audits remain required, and the coupled source risk must
not be replaced by the target's independent-block loss.

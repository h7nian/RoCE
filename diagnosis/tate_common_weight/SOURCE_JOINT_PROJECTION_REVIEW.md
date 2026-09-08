# Coupled source projection candidates

`source_joint_projection.R` extends the training-noise final-block candidates
with every initial-outcome and initial-density-ratio block. It follows the
backward adjoint equations and retains separate target/source calibration
contributions and actual source-piece fractions. Active M_fit truncation
derivatives are included; the target gradient summary remains untruncated,
matching the native calibration graph.

For each initial block, explicit per-observation coupling rows reconstruct
the corresponding Jacobian-transpose product. Their two-site variance scale
defines candidate L1 penalties through the shared noise-penalty helper. This
is a plug-in regularization proposal, not a concentration theorem: dependence
of the fitted downstream coefficients and Hessian fluctuations remains part
of the full statistical analysis.

On both first-validation-split systems, scales 0.5,1,2 each solved all eight
blocks. Maximum coupling reconstruction error is 1.25e-16 and maximum excess
of adjoint residual beyond its L1 band is 5.26e-8. Full coefficient L1 norms
are 4.928033/0.04598696/0 (control) and 5.761029/0.1008486/0 (treated).
No evaluation truth, coverage or validation risk was used to select among
these candidates. Validation views are explicitly rejected as fitting input.

This closes the omission of initial blocks in the candidate constructor; it
does not close source hyperparameter selection or inference. Selection must
respect the conditional coupled risks rather than a fictitious symmetric
joint quadratic. The other validation fits/audits and complete target/source
weight graph remain required before production adoption.

The constructor now accepts a separate `initial_penalty_scale`, defaulting
to the final scale for backward compatibility. Both arms' default scale-0.5
candidate objects are identical to the archived first-split objects. Changing
only the initial scale leaves final coefficients identical; explicit NULL
sets all initial coefficients to zero and labels them as null candidates,
not solved adjoint blocks. These separated-stage checks passed on both arms.
The fixed candidate selection procedure is recorded in
`SOURCE_STAGEWISE_SELECTION_PROTOCOL.md`; it has not yet been executed over
all validation folds.

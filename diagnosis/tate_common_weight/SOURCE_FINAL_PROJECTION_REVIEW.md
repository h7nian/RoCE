# Final source calibration projection: training noise scale

`source_final_projection.R` exposes the two final calibration blocks, not
the complete nuisance correction. It constructs separate target/source rows
of the final score derivative and checks their summed means against the
analytic gradient. Alpha has both target and source terms; gamma has only a
source term, including the active inference-truncation derivative.

For a block of dimension d, its coordinate penalties are

    c * sqrt(log(2*d) * (Var_target(row_j)/n_target
                         + Var_source(row_j)/n_source)).

This preserves separate site denominators. The constants {0.5,1,2} multiply
a coefficient penalty, not an SE. This training-only variance scale is a
candidate regularization rule, not a proven simultaneous concentration bound
under the full fitted/fold-balanced design. It must not be described as a
validated source projection policy.

Eight unit assertions pass, including unequal site sizes, site swapping,
location invariance, scale multiplication and malformed-input rejection.
On both audited first-split arm systems all 12 block/candidate solves
completed; maximum original-system KKT residual is 3.92e-8. Scale 2 produces
zero final coefficients in this check. This is not a reason to select or
reject that scale: no validation loss or TATE truth was used here.

The backward initial-fit corrections remain mandatory. A list of final-block
coefficients returned by this helper is deliberately not exposed as a full
TATE score or influence function. The remaining validation fits, source
conditional-risk policy and complete covariance analysis are still required.

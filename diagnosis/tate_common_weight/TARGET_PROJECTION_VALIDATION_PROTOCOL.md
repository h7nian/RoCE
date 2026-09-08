# Training-only target projection validation: next candidate protocol

This protocol is fixed before evaluating its candidate choices. It is still
method development informed by earlier seed10013 diagnostics; that seed's
outer fold must not be relabeled as independent validation. It introduces no
bootstrap or SE multiplier. No production estimator is selected yet.

## Fixed procedure for the first check

Use seed10013, outer fold1 and the captured target CV states. For each inner
validation fold k2=2,3,4,5, use the PS/treated-OR/control-OR models keyed by
(k1=1,k2), trained without either fold. Fit projection coefficients using only
the same out-of-two-fold training observations. Evaluate the coefficient risk
on k2, not on outer fold1. Do not use models trained on k2 to score that fold.

The three target moment Jacobian blocks are negative curvature matrices.
Writing H=-J and D for the target TATE score derivative, the unpenalized
coefficient equation is H a=-D. The validation loss is consequently

    0.5 a' H_validation a + D_validation' a.

This is a moment-projection risk, not squared error against the true TATE,
coverage, or interval width. Average validation losses using actual validation
fold sizes. A null projection a=0 is an explicit candidate, with loss zero;
the score need not receive a noisy correction when none is supported.

For nonnull candidates use the fixed constants c in {0.5,1,2}, with coordinate
penalties

    lambda_j = c * training_SD(individual score derivative_j)
                 * sqrt(log(2*d) / n_training).

Here d is the total target nuisance dimension (603 for this gate). These
penalties regularize the coefficient equations, not a reported SE. Unlike
fixed fractions of max|D|, they have a vanishing sample-size scale when the
coordinate variances remain controlled. That is rate-compatible motivation,
not a proof of the needed high-dimensional nuisance/projection bounds.
Penalties and all nuisance states for each split use training data only.

Choose the minimum validation-risk candidate; numerical ties prefer the null
projection, then the larger c. No one-standard-error coverage rule is used.
A nonnull candidate is ineligible if any required split fails its numerical
gate; retain every failure and do not average only successful folds. After
selection, refit the coefficient equations on the complete outer-training
set using its saved outer-complement nuisance models and the same c rule.

Outer fold1 is inspected only after this selection, reporting all candidate
risks, selected coefficients, score distributions and any large shifts. It
is not used to revise the constant grid or rank candidates. All such outputs
remain exploratory; an independently frozen validation experiment is required
before final adoption. The source ten-block graph, common-weight rule,
analytic variance and full simulation/RHC scope are not replaced by this
target-only coefficient-policy check.

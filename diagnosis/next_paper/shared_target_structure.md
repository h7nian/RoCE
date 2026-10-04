# A route from the reference model to calibrated source scores

Status: population derivation and implemented Gaussian reference, 2026-09-21.
This is not yet a growing-K theorem for the fitted RoCE procedure.

## A possible common target component

Consider one arm, a shared outcome basis, conditional-mean transportability
for the informative sources, unique population loss minimizers, and inactive
population truncation. The source outcome calibration loss uses the initial
merged weight; the checked implementation passes that plug-in through
`fit_unified_outcome` with the score-derivative weight convention.

If the outcome model is correct, the true outcome coefficient solves every
source's weighted outcome score, even if the initial weighting model is
misspecified. Uniqueness then gives a common final outcome function.

If the outcome model is misspecified but the initial merged weight is true,
the source population loss becomes, up to an irrelevant normalization,

    E_source[I(A=a) q_true {G(phi'alpha) - Y phi'alpha}]
      = E_target[G(phi'alpha) - m_target,a(X) phi'alpha].

This is the same target-population GLM projection for every informative
source. Thus common outcome working models and these correctness branches
can give a common source-side outcome limit even under outcome
misspecification. This argument does not apply merely because marginal
candidate biases cancel, nor does it cover arbitrary both-wrong fits or
persistently active truncation.

For source-assisted TATEs with both arms satisfying this structure, a
candidate first-order representation is

    Z_j - theta = U_target + V_source,j + remainder_j,

where U_target is the average of the common centered prediction contrast
over target observations, and V_source,j is the centered local residual
correction. Private source terms are independent across sites at the leading
population-score level. Their variances differ with the source law and its
sample size. The target anchor can be correlated with U_target, so the full
anchor/source covariance is not an independence matrix.

The next proof must establish the common population-limit argument under the
actual one/two-round loss definitions, and uniform negligibility of the
remainders across K. Current cross-fitted source predictions need not be
exactly equal in finite samples. A common-component approximation must not
be imposed on raw fitted outputs without bounding that difference.

## Known shared-component covariance with unequal private variances

Suppose the valid candidate errors are

    Z_j - theta = sqrt(h) U + sqrt(d_j) V_j,

with independent standard normals U,V_j, positive d_j, and marginal standard
errors sqrt(h+d_j). Invalid candidates may have arbitrary means; calibration
only needs to be correct on a subset of q valid candidates.

For a candidate interval using critical value c, conditional on U=u the
valid candidate's coverage probability is

    p_j(c,u) = Phi((c*sqrt(h+d_j)-sqrt(h)*u)/sqrt(d_j))
                - Phi((-c*sqrt(h+d_j)-sqrt(h)*u)/sqrt(d_j)).

At each u, choose the q smallest of the K computed probabilities. The
Poisson-binomial probability of at least r successes on this subset is no
larger than that for any fixed valid subset of size q, by coordinatewise
monotonicity. Integrating this conditional lower bound against phi(u) gives
a conservative lower bound for the r-vote event without enumerating all
subsets. The minimizing subset is allowed to vary with u only to obtain a
conservative bound; this is not a source-selection rule applied to the data.

The reference implements this one-dimensional calibration using a sorted
probability step and Poisson-binomial dynamic program. Tests verify reduction
to equicorrelated calibration when all d_j agree and to the independent
beta-order-statistic formula when h=0. The variance-region extension has a
separate proof and prototype in variance_uncertainty.md.

## A target-local estimator of the common projection

For arm a, consider a target weighted outcome fit minimizing

    E_target[I(A=a)/pi_t,a(X) * {G(phi'alpha)-Y*phi'alpha}].

When the target propensity is correct, this equals the same target-population
projection loss above, even if the outcome model is misspecified. When the
outcome model is correct, any positive limiting propensity weights preserve
the true outcome solution, subject to uniqueness. Thus this additional
target-local projection can have the common source outcome limit in either
correctness branch, provided target-PS correctness is included in the
weight-correct branch. Ordinary unweighted target outcome fitting need not
have that projection under outcome misspecification.

A cross-fitted prediction contrast from these two target-local projections
could estimate the common target variance h. Each source can estimate its
private residual variance from its own calibrated fits and local summaries.
This avoids estimating h by averaging outcome coefficients from potentially
invalid sources. It does not yet prove uniform consistency or covariance
inference: those rates and sample-splitting details remain necessary. The
additional target-local fits would not themselves require another source
communication round. Current production code is unchanged.

## Estimated covariance and growing-K work still required

- Estimate the common projection/target variance without presuming a known
  valid-source set. A target IPW outcome projection may provide the same
  population limit under the shared correctness branches, but this needs a
  separate sparse consistency argument.
- Estimate private residual variances from local summaries, preserving the
  source/target sample-size scaling and cross-arm covariance.
- Control the effect of estimating h,d_j on the critical value. It is unsafe
  to assume that a vote critical value is always monotone in correlation:
  the maximum-vote and partial-vote cases can behave differently.
- Use the independence of leading private source sums, condition on the
  target component, and establish uniform approximation errors. A generic
  full-vector total-variation bound may impose unnecessarily strong K/n
  restrictions; a vote-specific argument may be sharper.
- Distinguish the number of informative sources from a statistical set of
  accidentally unbiased candidates. State which set the new guarantee uses.
- Extend the scalar construction to arm-specific validity and TATE
  projection. The scalar prototype currently gives up useful control-arm
  information from a source whose treated arm is incompatible.

These steps are required progress toward the full objective; the existing
known-Gaussian experiment does not establish them. Compare explicitly with
Guo's two JRSSB papers and the existing RIFL construction before making any
novelty or efficiency claim.

# Uniform calibrated-score remainder: a sufficient modular route

Status: a conditional lemma for the full-method proof. The rate and population
orthogonality assumptions below still require verification for the actual
one/two-round calibrated fits. This is not an assertion that lambda-min CV
automatically supplies the required rates.

## Candidate, population limits and folds

For one arm and valid source j, write the candidate score as the sum of a
target prediction average and a private source correction average:

    Z_j = P_target m_hat_j + P_source,j [I(A=a) q_hat_j (Y-m_hat_j)].

All nuisance functions used on an evaluation fold exclude that fold at every
site. There are F fixed folds, with each site's fold sizes comparable to its
sample size. Assume n_target and n_j are comparable to n_min; unequal site
growth needs the corresponding separate rates instead of this simplification.
Each valid source has fixed population nuisance limits m_bar_j,q_bar_j that
identify the desired target mean. For the common-factor reduction additionally
require m_bar_j=m_bar across valid sources, as in shared_target_structure.md.

Let zeta denote finite-dimensional working-model coefficients. On an event
whose probability tends to one uniformly over all valid sources and folds,
assume:

1. Prediction-score differences have L2 norm at most C*r_n in each evaluated
   population and a uniformly bounded envelope.
2. The fitted coefficient error is in a region where the population score
   drift satisfies

       |Psi_j(zeta_hat_j)-Psi_j(zeta_bar_j)| <= C*r_n^2.

3. The ideal scores are bounded, have uniformly controlled standardized third
   moments and nondegenerate variance on the relevant scale.

Condition 2 follows from zero population gradients along the working-model
parameters plus a uniform quadratic Taylor bound and coefficient rates in
the corresponding Hessian norm. Prediction consistency alone is insufficient
under misspecification. The bounded-envelope assumption can instead be
replaced by suitable conditional Bernstein-type tail conditions, but that
extension is not supplied here.

## Remainder bound

Condition on one fold's training data. Its held-out target/source observations
are independent of the fitted nuisance functions. Bernstein's inequality and
a union bound over K sources and F folds give, with ell_n of order log(KF/delta),

    max_j |R_j| <= C * [r_n^2 + r_n*sqrt(ell_n/n_min) + ell_n/n_min]

outside the nuisance-rate failure event and an additional event of probability
at most delta. This applies to a fixed population valid set; independence of
the fitted models across sources or across folds is not needed for the union
bound. The same argument handles the two arms together, changing constants.
Combining F fold contributions preserves the order when F is fixed.

The ideal source-assisted TATE contribution is then a common centered target
prediction contrast plus an independent centered source residual contrast.
The normalized uniform remainder is negligible if

    sqrt(n_min)*r_n^2 -> 0,
    r_n*sqrt(ell_n) -> 0,
    ell_n/sqrt(n_min) -> 0.

For a demonstrable sparse rate r_n=O(sqrt(s*log(pKF/delta)/n_min)), a sufficient
combined condition is

    s*log(pKF/delta) = o(sqrt(n_min)),

with s>=1 and the same fixed-F/comparable-size assumptions. This is a
sufficient regime, not a necessary or optimal condition.

## Why calibrated limits can be orthogonal

For a canonical outcome link m_alpha and log-linear merged weight q_beta,
the population coefficient gradients of the candidate mean are

    d_alpha Psi = E_target[dot_m_alpha * phi]
                 - E_source[I(A=a) q_beta dot_m_alpha * phi],
    d_beta Psi  = E_source[I(A=a) q_beta * phi * (Y-m_alpha)].

Here q_beta=exp(phi'beta). The implementation uses the opposite-sign tilt
coefficient gamma, so beta=-gamma; its gradient has the opposite sign and
the same zero-moment condition.

If both models are correct, both gradients vanish by transport and conditional
mean identities. If the outcome model is correct but the weight model is
misspecified, the second gradient vanishes; the first requires the weight
calibration moment with the correct limiting outcome derivative. If the
weight model is correct but the outcome model is misspecified, the first
gradient vanishes; the second requires the weighted outcome calibration
moment. These are population statements, not finite-sample exact equalities.

This identifies specific checks for the fitted algorithm: correct arm coding,
feature bases and normalizations; the derivative convention; initial fits
converging to the required branch; inactive population clipping; and final
losses having the claimed population minimizers. One-round target initial
outcome fits and two-round source fits cannot be interchanged without checking
these conditions. The existing shared-target argument is conditional on them.

## Combining with confidence-set inference

For the shared-target Gaussian approximation in growing_k_score_reduction.md,
bounded standardized third moments give the additional sufficient restriction
q/sqrt(n_min) -> 0 through the elementary sum of scalar approximation errors.
That restriction is conservative. The arbitrary-dependence marginal reference
has a different error bound and should be analyzed separately.

The variance-region inference also needs an actual simultaneous region from
estimated scores, with an explicit allowance for nuisance estimation error.
Sample-variance concentration for oracle scores alone is not enough. Nor is
an unproved O_p remainder an implementable radius inflation: a deterministic
high-probability bound, a valid estimate, or a justified asymptotically
negligible perturbation argument is required.

Thus the remaining proof tasks are concrete: establish uniform calibrated-fit
rates and gradients for each correctness branch, verify the common target
limit, construct fitted-score variance regions, and justify the resulting
confidence-set approximation. These requirements preserve the original goal
of a high-dimensional, weak-deviation, growing-K method rather than treating
the Gaussian reference experiment as the completed extension.

# Sufficient conditions for two-layer weight-learning moments

This note completes the logical route for the proposed outer-fit reuse rule.
It is a sufficient-condition argument, not a claim that the current manuscript
has already proved every condition for its cross-validated penalized fits.
The current study fixes K and F; n below denotes the common order of the
per-site training sample sizes.

## Means without an independent eta-validation fold

Fix outer fold k. Let M_T(nu) be a source-assisted mean, target mean, or their
difference, evaluated on the outer training data with nuisance parameter nu.
M_T combines the appropriate site-specific empirical means. Let M(nu) be its
population counterpart, nu_star the limiting nuisance target, and Delta the
error of the final outer-calibrated fit.

For a valid candidate, assume:

1. M(nu_star) is the intended mean/difference and its full nuisance gradient
   is zero at nu_star. This is orthogonality of the COMBINED estimating
   function. The target and source derivatives need not vanish separately.
2. The empirical gradient error at nu_star is O_p(sqrt(log(p)/n)) in sup norm,
   uniformly over the fixed collection of sites, arms and outer folds.
3. The final nuisance L1 error is O_p(s sqrt(log(p)/n)).
4. The empirical Taylor remainder after the linear term is
   O_p(s log(p)/n), for example through empirical prediction-error bounds and
   controlled weighted Hessians.
5. The limiting score has the moments needed for its empirical CLT.

The Taylor decomposition, applied on the same training observations, is

    M_T(nu_hat) - M(nu_star)
      = [M_T(nu_star) - M(nu_star)]
        + [gradient M_T(nu_star) - gradient M(nu_star)]' Delta
        + remainder.

Holder's inequality bounds the middle term by the product of the empirical
gradient sup norm and the nuisance L1 error. It does not require independence
of Delta and the training sample. Consequently,

    M_T(nu_hat) - M_T(nu_star) = O_p(s log(p)/n).

If s log(p) = o(sqrt(n)), this is o_p(n^(-1/2)). In particular, valid
source/anchor training discrepancies remain O_p(n^(-1/2)). All training
observations can contribute even though the final fitted nuisances reuse them.

The conditions must be checked for the actual calibrated score and limiting
parameters. Cross-fitting the INITIAL models inside the calibration loss does
not itself prove orthogonality at the FINAL nuisance limits. Cross-validation
and numerical fitting also need to deliver the stated rate events. Bounded
features/outcomes and a fixed weight radius help control empirical gradients,
but do not replace those estimation arguments.

For a separated biased candidate, mean consistency is sufficient for its
estimated discrepancy to converge to a nonzero bias. Root-n regularity of every
biased candidate is not required merely to force its Wald discrepancy to
diverge. This is distinct from the next-paper all-source Gaussian procedures,
which may require regularity of biased candidates too.

## Covariances

For each site let U_hat and U_star denote its fitted and limiting vector of
score contributions, with fixed dimension when K is fixed. Assume empirical
L2 consistency, ||U_hat-U_star||_(T,2)=o_p(1), and O_p(1) empirical second
moments of U_star. By Cauchy-Schwarz,

    |P_T(U_hat,j U_hat,l - U_star,j U_star,l)|
      <= ||U_hat,j-U_star,j||_(T,2) ||U_star,l||_(T,2)
         + ||U_hat,l-U_star,l||_(T,2) ||U_star,j||_(T,2)
         + ||U_hat,j-U_star,j||_(T,2) ||U_hat,l-U_star,l||_(T,2)
      = o_p(1).

Together with the limiting-score LLN and consistency of centering means, this
gives consistent within-site covariance matrices. Combine site components with
the actual sample-size weights; preserve cross-arm covariance. The same
argument applies to the implemented fixed-F block-centering convention because
the limiting block means agree and every block size diverges.

For a biased score this establishes covariance of its limiting score values,
not automatically the full sampling variance of its fitted estimator. The
fixed-K oracle aggregation argument ultimately uses the valid-source subspace.

## Eta and final outer inference

With consistent covariance coefficients, a stable unique oracle minimizer,
separated biases, and the manuscript's growing-cutoff regime, the fixed-K
selection/argmin argument yields eta_hat^(-k) -> eta_star and removes biased
coordinates with probability tending to one. On that event, each valid
outer-fold difference D_(j,k) is O_p(n^(-1/2)), so

    sum_j [eta_hat_j^(-k)-eta_star_j] D_(j,k) = o_p(n^(-1/2)).

No extra independent eta-validation fold is required for this product argument.
It is enough to prove the preceding training-moment bounds for the two-layer
procedure. The fixed cutoff2 used in simulation does not fulfill the
growing-cutoff assumption, so this argument does not certify its finite-sample
coverage or imply that eta sensitivity can always be omitted there.

A simple example explains the cutoff distinction. Suppose the one-source
rescaled variance objective tends to H eta^2 - 2 b eta + C, with H,b>0,
and the source is valid. Its Wald statistic can converge to |Z|, where
Z is standard normal. With a fixed lambda, minimizing this quadratic plus
(lambda |Z|-1)_+ |eta| gives the random limit

    eta_limit = max(0, [b - (lambda |Z|-1)_+/2]/H).

There is positive probability of both an unpenalized and a shrunk weight.
Thus constant-cutoff weight consistency is not automatic even when every
candidate is valid. This does NOT by itself prove invalid inference: for an
independent single holdout, conditional inference given the learned weights
may work without deterministic weight convergence. But the oracle replacement
proof above cannot be invoked for that setting, and the full cross-fold
estimator needs its own dependence and studentization argument.

## What changes when K grows

Uniform rather than pointwise score and covariance control becomes necessary.
For illustration, if max_j |D_j| = O_p(sqrt(log(K)/n)), the weight replacement
error is bounded by

    ||eta_hat-eta_star||_1 O_p(sqrt(log(K)/n)).

The required bound must be relative to the actual oracle standard error s_N,
not an automatically inherited sqrt(N_all) scale. If s_N is of order n^(-1/2),
a sufficient condition is ||eta_hat-eta_star||_1 sqrt(log(K))=o_p(1). If the
oracle gains a full K-fold variance reduction, so s_N is of order (nK)^(-1/2),
the corresponding sufficient condition strengthens to
||eta_hat-eta_star||_1 sqrt(K log(K))=o_p(1). These illustrative maximum-score
rates themselves need appropriate uniform tail/remainder assumptions.

With n observations per site, site fractions equal 1/(K+1); they vanish when
K grows. The present paper's nonvanishing-site-fraction regime therefore cannot
be carried over unchanged. Three-layer splitting can help organize independent
validation arguments but does not establish these rates or resolve weak-bias
selection by itself.

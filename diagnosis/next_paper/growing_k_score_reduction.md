# A modular coverage bound for independent site score sums

Status: a reduction under explicit score and approximation assumptions. This
does not yet verify those assumptions for every fitted RoCE component.

Fix a population subset V of q informative sources. Suppose their candidate
estimates have representations

    Z_j - theta = U + V_j + R_j,  j in V,

where U is a centered target score average, the V_j are centered private
source score averages, and U,V_1,... are mutually independent. Dependence of
the target anchor on U is permitted. Let h=Var(U), d_j=Var(V_j)>0 and
sigma_j=sqrt(h+d_j). The remainders R_j may be dependent.

## Gaussian approximation without a full K-dimensional CLT

For h>0 assume the Kolmogorov distance of U/sqrt(h) to N(0,1) is at most e_U.
If h=0 and U=0 almost surely, set e_U=0. Suppose the corresponding distance
for V_j/sqrt(d_j) is at most e_j. For any fixed critical value c, conditional
coverage of a valid source differs from its Gaussian counterpart by at most
2e_j, uniformly in the value of U.

Conditional on U, source coverage indicators are independent Bernoulli
variables. Coupling these Bernoulli variables bounds the change in the
probability of at least r successes by 2*sum_{j in V} e_j.

The Gaussian conditional r-vote probability is an even, nonincreasing
function of |U|: every symmetric-normal interval probability has that
property, and a Poisson-binomial tail is increasing in each success
probability. This function has total variation at most2. Replacing U by its
Gaussian counterpart therefore changes its expectation by at most2e_U,
using integration by parts against the two distribution functions.

Thus the actual ideal-score r-vote probability is at least its Gaussian
counterpart minus

    2e_U + 2*sum_{j in V} e_j.

This comparison concerns a fixed truly informative subset V. The algorithm
need not identify it; the calibration already lower-bounds every such subset.

## Remainder and variance events

Let A_var be a simultaneous variance-region event as in
variance_uncertainty.md, with failure probability delta_var. Let A_rem
satisfy, for the anchor and every informative source,

    |R_0| <= a_0,n * sigma_0,
    |R_j| <= a_n * sigma_j,

with failure probability at most delta_rem. The tolerances are part of the
procedure and must dominate genuine uniform remainder bounds. Their choice
is not justified by observing good simulation coverage.

Inflate the standardized source radius by a_n and the anchor radius by
a_0,n, using the variance upper bounds. On A_var intersect A_rem, the
reported intervals contain the corresponding ideal-score intervals formed
with the true variances and fixed population critical values. If the ideal
anchor's normal approximation error is e_0, the confidence-set coverage is
at least

    1 - alpha_A - alpha_S - delta_var - delta_rem
      - e_0 - 2e_U - 2*sum_{j in V} e_j.

Hulling the confidence set and the declared empty-set fallback do not lower
this coverage bound. Independence of the covariance estimates, fitted
remainders and mean estimates is not assumed in this event argument.

## One sufficient growth regime, and its limits

For independent within-site score observations with uniformly bounded
standardized third absolute moments, scalar Berry–Esseen bounds give
e_j=O(n_j^{-1/2}), e_U=O(n_t^{-1/2}), and e_0=O(n_t^{-1/2}). A conservative
sufficient source-count condition is therefore

    q_n / sqrt(n_min) -> 0,

or K_n=o(sqrt(n_min)) when q_n is proportional to K_n. This bound is not
claimed to be optimal. It avoids silently treating a fixed-K theorem as a
growing-K theorem, but faster K growth would require a sharper vote-specific
comparison or stronger distributional assumptions.

To obtain shrinking intervals comparable to the known-score reference also
require a_n,a_0,n -> 0, vanishing event failure probabilities, and suitable
variance-region contraction. None of those follows solely from this lemma.
Uniform high-dimensional nuisance rates, source/target sample-size scaling,
and validity of variance bounds from fitted scores remain required work.
The current lambda-min cross-validation policy is not asserted to select
theoretical regularization rates automatically.

The arbitrary-dependence marginal reference in theory.md has a different
error bound, proportional to q_n/(q_n-r_n+1) times marginal approximation
errors. It may allow weaker growth restrictions at the cost of wider
intervals. The two bounds should be compared rather than conflated.

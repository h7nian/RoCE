# A computable subset-variance bound with private arm correlation

This is a proposed extension of the existing shared-prediction Gaussian
reference. It addresses a concrete issue in the fitted-score packets: their
private cross-arm covariance is retained and generally nonzero. It does not
resolve conditional nuisance bias or estimated-covariance inference.

Suppose source candidates share the same two target prediction components,
and independent source-private noise blocks have covariance

    V_j = [[v_j1, c_j], [c_j, v_j0]],
    rho_j = c_j / sqrt(v_j1 v_j0).

For positive marginal variances and a positive-semidefinite V_j,

    V_j <= Vbar_j = (1+abs(rho_j)) diag(v_j1,v_j0)

in Loewner order. Indeed the difference has nonnegative diagonal and zero
determinant. This inflation removes the private cross-arm correlation while
dominating its entire covariance block, irrespective of its sign.

Keep the original target covariance and replace only private blocks to form
Sigma_bar >= Sigma. Every valid-coordinate subset retains this order. Its
minimum variance among unbiased linear estimators of the TATE therefore
satisfies

    v_V(Sigma) <= v_V(Sigma_bar).

The inflated covariance has the diagonal-private structure of the existing
shared_prediction_exact calculation. Its worst subset at prescribed arm
counts can be found by choosing the largest inflated private variances in
each arm. Thus its maximum subset variance is an upper bound for the original
problem without enumerating all validity patterns.

The required bias coefficient is

    C_upper = min(max_V v_V(Sigma_bar), v_known_reference(Sigma))
              - v_all(Sigma).

Every allowed subset contains the known-valid reference, so its ORIGINAL
variance provides an additional upper bound. Taking the minimum avoids
weakening the existing known_valid_upper bound when inflation also affects
source arms already known to be valid.

The subtraction must use the ORIGINAL full-pool variance. Subtracting the
inflated full-pool variance instead can invalidate the bound. The full-pool
GLS weights, its variance, the reference estimator and the noncentrality
statistic should likewise retain the original covariance; inflation is only
used to upper-bound the worst valid-subset variance.

When all private correlations are zero this reduces to the existing exact
shared-prediction calculation. Cases in which every coordinate is guaranteed
valid should retain their direct GLS interval without an unnecessary bias
bound. Degenerate covariance blocks require explicit handling or rejection;
silently applying a numerical inverse is not a validity argument.

With bounded comparable private variances and q_a proportional to K, the
inflated subset variance and full-pool variance approach the same shared
target floor. Their difference can remain O(1/(nK)), preserving the relevant
Gaussian fourth-root adaptation scale up to constants. A full statement needs
the target/anchor covariance and nondegeneracy conditions; this paragraph is
not a fitted-method rate theorem.

The bound_method argument now includes shared_prediction_bound. Check the bound
against exhaustive subset enumeration for small K, positive and negative
private correlations, different valid counts and both arm patterns. Compare
the zero-correlation case with the existing reference and check actual fitted
population packets. The initial checks cover those cases; the extension is
still a Gaussian-reference component, not a validated fitted-method interval.

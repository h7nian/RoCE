# Common predictions and the population weight problem as K grows

This is an algebraic observation about the oracle score covariance, not a
completed selection or inference theorem. It uses the population conditions
in common_outcome_score_orthogonality.md and equal per-site sample size n.

For a fixed arm, each source candidate has the form C+R_j, with the SAME
target component C and an independent source residual mean R_j. Let T be
the target anchor. If S is a set of q valid sources, define

    s = sum_{j in S} eta_j,
    v_j = Var(R_j),
    P = sum_{j in S} 1/v_j.

The aggregate is (1-s)T+sC+sum eta_j R_j. At a fixed s, its private variance
is minimized by eta_j=s/(P v_j), giving s^2/P. Consequently the oracle
variance problem reduces to one scalar sum s and

    V(s) = (1-s)^2 Var(T) + 2s(1-s) Cov(T,C)
           + s^2 Var(C) + s^2/P.

At an unconstrained valid-source optimum, the common part of the derivative
with respect to any source coefficient equals -2s/P. An invalid source with
coefficient zero has this same smooth derivative: its private residual is
independent of the currently used valid-source residuals. If |s| is bounded,
q>=cK, and v_j<=C_v/n, then on the existing N_all=(K+1)n objective scale,

    N_all |partial_eta_bad V| <= 2 (K+1) |s| C_v/q = O(1).

This removes one possible K-dependent contribution from source-specific
target prediction functions. In the original construction, an invalid
source's target prediction can differ from every valid prediction. Its
covariance with the shared part of the oracle estimator need not decrease
with K, so multiplying that derivative by N_all can instead give an O(K)
term. This is a worst-case distinction, not a claim about every current DGP.

## Two-arm extension and its scope

At the relevant population nuisance limits, each valid-arm residual has mean
zero. This follows from the correct OR, or from correct transport weights and
the intercept equation of the full-target IPW OR projection. Since the two
arm residuals from one person cannot both be nonzero, their covariance is
minus the product of their means. Thus the private cross-arm covariance is
zero whenever at least one of those arms is valid, at these population limits.

The valid oracle problem then reduces to two sums s1/s0, a 2-by-2 target
covariance calculation and arm-specific precision sums P1/P0. The same
O(K/q_a) derivative bound applies to an excluded coordinate of arm a, if the
two-sum optimum is uniformly bounded and q_a is a positive fraction of K.
An invalid partner arm has oracle weight zero; a valid partner has zero
population private cross-covariance with it.

Finite fitted-score residual means are generally NOT zero conditional on
training. Their cross-arm covariance must still be estimated and retained.
The implementation checks explicitly retain it. Proving that replacement by
the population structure is negligible requires uniform nuisance and covariance
rates at the scale of the desired interval, not just fixed-K consistency.

## Relation to separated and weak departures

For an ideal uniform Gaussian-type Wald approximation, making valid-source
penalties vanish uniformly suggests a cutoff c_nK much larger than
sqrt(log K). Making penalties diverge for all separated invalid arms suggests
c_nK much smaller than sqrt(n) times the minimum absolute bias. These two
requirements need a separation window; they cannot both hold for biases of
order n^-1/2 as K increases.

The common-prediction structure can simplify the separated-source weight
problem, but it does not solve weak-deviation inference. Confidence intervals
must still allow uncertain source validity and valid-candidate nuisance bias.
The current K4/6 experiments vary fixed K, and their rho=.5 departures are
finite-sample boundaries, not an asymptotic local-bias sequence.

The same factor structure may permit a small-dimensional or blockwise weight
solver, but no computational complexity or finite-sample improvement has been
validated for a dedicated implementation yet. No current-paper default has
been changed on the basis of this observation.

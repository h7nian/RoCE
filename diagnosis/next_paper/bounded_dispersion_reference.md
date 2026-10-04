# Replicated bounded observations: a finite-sample reference

This construction uses independent observations within and across sites,
known bounds on their ranges, and a prespecified minimum number of valid
source means. It needs neither Gaussian noise nor known outcome variances.
The interfaces cover scalar means and a two-arm contrast, with a joint error
allocation that does not assume independence between arms. This is not yet a
high-dimensional fitted-score theorem or a finite-sample efficiency claim.

Let source j contain n>=2 observations Z_ji of mean mu_j and variance v_j,
with range width at most R. Sites may have different distributions. At least
q source means equal the target theta. Put

    B = sum_j(mu_j-mubar)^2/(K-1),
    Bhat = mean_j[(S_j^2-T_j)/(n(n-1))]
           - [(sum_j Zbar_j)^2-sum_j Zbar_j^2]/[K(K-1)],

where S_j=sum_i Z_ji and T_j=sum_i Z_ji^2. Only sums and squared sums need
to be transmitted. E Bhat=B; the within-site term removes measurement-noise
variance rather than interpreting all observed dispersion as heterogeneity.

## Exact variance bound

Write Bhat=Z'AZ. Matrix A has zero diagonal, off-diagonal entries
1/[K n(n-1)] within a site, and -1/[K(K-1)n^2] between sites. Its row sums
are zero. For the repeated mean vector mu,

    ||A||_F^2 = 1/[K n(n-1)] + 1/[K(K-1)n^2],
    ||A||_op = 1/[(K-1)n],
    ||A mu||_2^2 = B/[(K-1)n].

Independence and the zero diagonal imply that the centered quadratic and
linear terms have zero covariance, without any symmetry assumption. Hence

    Var(Bhat) = 4 sum_i (A mu)_i^2 v_i
                 + 2 sum_(i,j) A_ij^2 v_i v_j
               <= a B+b,
    a = 4 v_max/[(K-1)n],
    b = 2 v_max^2 ||A||_F^2,
    v_max = R^2/4.

Cantelli's inequality gives, with probability at least 1-delta,
B-Bhat <= sqrt(c(a B+b)), where c=(1-delta)/delta. A computable upper bound is

    U = max{0, Bhat+c a/2
                + sqrt(max{0,c a Bhat+(c a)^2/4+c b})}.

If B>Bhat, squaring the concentration inequality gives its upper quadratic
root; if B<=Bhat, the same expression is at least Bhat. A negative
discriminant has no feasible positive B on the concentration event.

## Coverage and length

The validity-count inequality is

    |mubar-theta| <= sqrt{(K-q)(K-1)/(qK) * B}.

Add this bias allowance with B replaced by U to a Hoeffding interval for
the pooled source mean. Intersect it with a separately valid target interval,
using the shorter interval when the intersection is empty. The sum of the
dispersion, source-mean and target error budgets is alpha. Coverage is at
least 1-alpha for every set of source distributions satisfying the assumptions.
When q=0, use target-only inference; when q=K, pool the valid observations
directly. All constants in this reference are explicit.

At the all-valid law, B=0 and Var(Bhat)<=b. The upper root satisfies
U <= 2(Bhat)_+ + 2ca + sqrt(cb). Thus

    E U = O(R^2/(n sqrt(K)) + R^2/(nK)),
    E sqrt(U) = O(R/(sqrt(n) K^(1/4))).

For fixed q/K bounded away from zero and one, the resulting expected interval
length is O(R/(sqrt(n)K^(1/4))) at the all-valid law, and never exceeds the
target reference length O(R/sqrt(n)). This matches the lower-bound order
in the randomized-binary submodel, up to constants. The elementary Cantelli
calibration can be quite conservative at practical K, so this reference must
not be presented as already competitive with Gaussian or fitted RoCE intervals.

For the observed randomized binary model, Z=A Y/e is bounded in [0,1/e]
and has mean E[Y(1)]. Control scores have the analogous form. Source validity
is equality of the relevant potential-outcome means. Within-site replication
is essential to estimate the dispersion without knowing the noise variances.
Actual fitted transport scores may share target data and nuisance fits;
their dependence and approximation errors must be addressed before applying
this independent-score reference to the full method.

Checks verify the summary/quadratic identities, eigenvalues,
exact Bernoulli mean and variance calculations, one-sided coverage inversion,
and interval special cases. No main estimator or default DGP was changed.

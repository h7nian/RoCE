# Estimated variance: a sufficient coverage construction

Status: mathematical construction under a Gaussian shared-target model.
The Gaussian reference is implemented in variance_region_calibration.R.
Simultaneous variance bounds and the full fitted-score approximation have
not yet been established for the high-dimensional RoCE pipeline.

## Why simply substituting estimates is insufficient

The known-variance prototype chooses a critical value from a fixed covariance
model. A data-dependent covariance estimate changes that calibration. Even a
collection of intervals each having correct coverage at a fixed covariance
does not automatically justify selecting among them using the same data.
Neither general correlation monotonicity nor independence between fitted
means and variance estimates should be presumed.

## Simultaneous variance region

Let good source errors have variances h+d_j, with a common Gaussian target
component of variance h and independent Gaussian private components d_j.
Suppose an event A, with probability at least1-delta, provides

    h in [h_L,h_U], d_j in [d_L,j,d_U,j] for every valid source,
    Var(Z0) <= v_U,0.

The intervals may be data-dependent. Coverage of the variance region is only
required for the valid sources, but a construction using all-source
simultaneous bounds can supply it without knowing their identities.
Invalid-source variance bounds need not be correct for this coverage argument.

For a valid source, its standardized common-factor loading belongs to

    a_j in [ell_j,u_j],
    ell_j=sqrt(h_L/(h_L+d_U,j)),
    u_j=sqrt(h_U/(h_U+d_L,j)).

Its standardized error can be written a_j U + sqrt(1-a_j^2) V_j.

## A computable conditional lower bound

For a fixed common noise value z, the absolute conditional mean is at most
u_j |z|, and its private standard deviation lies between
sqrt(1-u_j^2) and sqrt(1-ell_j^2). For a symmetric interval [-c,c], normal
coverage decreases with the absolute mean. At a fixed mean, its minimum
over a closed standard-deviation interval occurs at an endpoint: if the mean
is inside [-c,c], coverage decreases with the standard deviation; if outside,
the coverage function has at most one interior stationary point, a maximum.

Consequently the smaller of the two endpoint normal probabilities at mean
u_j|z| is a lower bound for the valid source's conditional coverage. This
relaxes the coupling between its mean and standard deviation and can be
conservative. Zero private-variance endpoints need their limiting treatment.

Use the q smallest of these lower probabilities in the Poisson-binomial
tail calculation, then integrate over z. The result lower-bounds the r-vote
CDF for any fixed valid subset of size q and every covariance inside the
region. Calibrating this lower bound at1-alpha_S therefore gives a common
standardized critical value c_hat that dominates the true fixed-subset
critical value on A. The universal marginal critical value is also valid;
taking the smaller of two sufficient critical values preserves this
domination. If the covariance envelope cannot improve the marginal bound,
the marginal bound is the declared fallback.

The distinction between standardized calibration and calibration using
arbitrarily enlarged, data-dependent physical radii matters. The latter
requires its own argument and must not be treated as automatically valid.

## Confidence statement

Use source radii c_hat*sqrt(h_U+d_U,j), which dominate the radii based on true
source standard errors on A. Use anchor radius
z_(1-alpha_A/2)*sqrt(v_U,0). Form the same quorum confidence set intersected
with the anchor interval. For alpha_A+alpha_S+delta=alpha, coverage is at
least1-alpha under the Gaussian model and the simultaneous-region guarantee.

To see why dependence on the data is allowed, fix any q truly valid indices.
The true quantile for their standardized vote statistic is deterministic
under the data-generating law. On A, the constructed radii dominate those
using this fixed quantile and the true variances. Failure is thus contained
in the union of A-complement, the fixed anchor tail event, and the fixed vote
tail event. Their probabilities are bounded by delta,alpha_A,alpha_S.
No conditioning argument asserting independent variance estimates is used.

## Prototype and full-method requirements

An unknown-variance Gaussian-score experiment can construct exact
chi-square intervals from independent normal score samples for the shared
and private variance components, with a simultaneous probability allocation.
This is still an oracle-score experiment. It does not validate variance
regions made from estimated nuisance scores.

The implementation checks the endpoint lower bound, degenerate limits,
known-variance reductions, probability allocation and scale equivariance.
The initial infinite-domain quadrature and then fixed finite pieces failed
on kinked conditional probabilities. The revised integrator bisects failed
pieces while preserving the total absolute-error budget. If that refinement
still fails, calibration retains the proved marginal rule and records the
numerical fallback reason; it never accepts a failed quadrature value.
The omitted standard-normal mass outside [-8,8] is less than 2e-15 and its
omission decreases the computed lower bound. Failed inputs and test versions
are retained in scratch implementation/next_paper/v3.

Before a full-method claim, establish valid simultaneous regions for actual
estimated scores and bound nuisance remainders uniformly as K and p increase.
Broader Monte Carlo calibration and larger-K computational checks also remain.
The small Gaussian pipeline checks are implementation checks, not precise
coverage estimates. No current-paper estimator is changed by this prototype.

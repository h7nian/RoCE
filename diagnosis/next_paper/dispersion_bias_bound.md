# A dispersion-based bias bound: next candidate after the quorum reference

Status: a derived Gaussian-reference route. The covariance identity, its
computable upper bound and the bias inequality passed280 random subset checks
(maximum identity error1.03e-13). It is not a novelty claim or a validated fitted
RoCE procedure. The failed efficiency of the interval-voting reference
motivates testing this alternative.

The known-Gaussian interval is now implemented in `dispersion_bias_intervals.R`.
Its59 checks cover noncentrality inversion, exact and conservative covariance
bounds, known valid anchors and two-arm variance. The24x1000 paired reference
study is complete. Against a baseline exploiting the same fully-valid-arm
assumption, the new anchored interval is still wider in the examined
treated-half cases. A general efficiency solution remains open.

## Scalar homoscedastic case

Let X_j=theta+b_j+epsilon_j, j=1,...,K, with independent N(0,d) errors and
at least q>=1 zero biases. Put gamma=q/K and m=K-q. Cauchy-Schwarz on the
nonzero bias coordinates gives

    bbar^2 <= m/q * (1/K) sum_j (b_j-bbar)^2.

For Q=sum_j(X_j-Xbar)^2/d, Gaussian projection gives a noncentral chi-square
law with K-1 degrees of freedom and

    lambda=sum_j(b_j-bbar)^2/d.

Let U(Q) be a one-sided 1-delta upper confidence bound on lambda obtained
by inverting this distribution. Then a valid reference interval is

    Xbar +/- [z_(1-alpha_noise/2)*sqrt(d/K)
              + sqrt{d*(K-q)/(K*q)*U(Q)}],

with alpha_noise+delta=alpha. A union bound proves coverage; independence of
the Gaussian mean and residual sum of squares also holds in this model.
If q=K, the bias bound is exactly zero and no dispersion budget is needed.

At the all-valid point, U(Q)=O_p(sqrt(K)) for fixed delta. For q/K bounded
away from0 and1 the bias allowance is O_p(sqrt(d)*K^(-1/4)). With d proportional
to1/n this suggests the n^(-1/2)K^(-1/4) local adaptation order discussed in
`majority_oracle_adaptation_bound.md`. This rate statement needs an explicit
expectation/tail argument before it can be called a matching expected-length
upper bound. It does not contradict the half-valid lower bound at ambiguous
two-group configurations, where dispersion and the interval remain large.

An additive common N(0,h) error changes the Gaussian mean variance to
h+d/K but cancels from Q and the bias bound. Thus a nonvanishing shared-target
variance still places a floor on the final interval length.

The pooled-source bias allowance need not shrink under fixed, strong
departures: genuine between-source bias dispersion then remains nonzero.
Consequently this interval alone is not an adequate final estimator/inference
solution over the whole model class. Intersecting it with a separately valid
anchor interval, with a joint error allocation, can retain target-rate
consistency while allowing improvement near the all-valid point. Oracle
efficiency under separated alternatives is a further requirement to study,
not something established by this dispersion bound.

## Known general covariance

Let X~N(theta*1+b,Sigma), Sigma positive definite. Define

    v = 1/(1' Sigma^(-1) 1),
    w = v*Sigma^(-1)1,
    P = Sigma^(-1) - Sigma^(-1)1*v*1'Sigma^(-1).

Then M=w'X has noise variance v, Q=X'PX has a noncentral chi-square law with
K-1 degrees of freedom and lambda=b'Pb, and their Gaussian noise components
are independent.

For a fixed zero-bias subset V, let v_V=1/(1_V' Sigma_VV^(-1) 1_V).
The exact quadratic bound is

    (w'b)^2 <= (v_V-v) * b'Pb,  whenever b_V=0.

One derivation constrains r=b-(w'b)1 by w'r=0 and r_V=-(w'b)1_V.
Minimizing r'Sigma^(-1)r under these linear constraints gives
(w'b)^2/(v_V-v). The q=K case has b=0 and is handled separately.

Consequently C_q=max_{|V|=q}(v_V-v) gives a valid bias allowance sqrt(C_q*U).
For a fixed support S=V-complement, the equivalent identity is

    w_S' (P_SS)^(-1) w_S = v_V-v.

Exact subset enumeration is feasible only for small K. A computationally
cheap upper bound follows from uniform weights on V:

    C_q <= max_i Sigma_ii/q + (q-1)/q * max_(i!=j) Sigma_ij - v.

Use zero when q=K. The bound can be loose for general heterogeneous target
loadings; it is exact for compound-symmetry covariance. For
Sigma=h*11'+diag(d_j), exact worst-subset computation selects the q largest
d_j and gives C_q=1/sum_(q largest d_j)(1/d_j)-1/sum_all(1/d_j).

## Different valid arm sets

Apply the bias-dispersion bound separately to the two source mean vectors,
using their prespecified q1 and q0. Their valid subsets need not overlap.
Use the FULL joint arm covariance for the noise variance of M1-M0. A valid
union-bound interval has radius

    z*sd(M1-M0) + sqrt(C1*U1) + sqrt(C0*U0),

with the two dispersion failure budgets plus the noise budget totaling alpha.
If an arm is assumed fully valid, its bias allowance is zero. Jointly
optimizing weights is a separate refinement; the first reference uses the
explicit marginal GLS means and their actual cross-arm covariance.

If a trusted anchor is included in an arm's candidate vector, its index is a
known valid coordinate. The worst-subset bound must maximize only over valid
subsets containing that coordinate. In particular C<=v_anchor-v is always
available. This extends the reference to zero guaranteed valid sources without
silently assuming a valid source exists; the anchor still supplies one valid
candidate. A growing-K efficiency bound for this constrained construction
depends on its covariance structure and must be derived separately.

## Assumptions that must not be hidden

Unlike the quorum proof, the noncentral chi-square pivot requires a joint
Gaussian approximation around the biased means of ALL participating sources,
plus adequate covariance control. It does not allow arbitrary invalid-source
error distributions. A valid-source-only Gaussian approximation is insufficient.

This matters for our one-round initialization: under an OR-correct but
weight-misspecified source with a shifted outcome, target-initialized OR may
not equal that source's true initial OR. A root-n expansion around its biased
limit cannot be assumed from the valid-source theorem. Two-round source-local
initialization, or a suitable honest evaluation construction, needs auditing
before claiming that this reference applies to every fitted candidate.

The population audit in `biased_source_initialization_audit.md` confirms a
nonzero derivative in shifted C3 sources with target initialization and a zero
derivative with source initialization. It also distinguishes fixed departures
from local sequences and the identity-link/correct-weight exceptions.

Estimated covariance, growing-K chi-square approximation, nuisance remainders,
and the error from using a common-target covariance structure remain open.
They must be proved or controlled explicitly. The current paper's estimator
and its running experiments are unchanged.

There is an additional rate issue for the quadratic statistic. If each
standardized candidate remainder is bounded by r, its projected Euclidean
norm can be of order sqrt(K)*r. Near the all-valid point Q is of order K,
so a direct perturbation bound for Q is of order K*r. Making this negligible
relative to the O(sqrt(K)) noncentrality allowance can require sqrt(K)*r=o(1),
which is stronger than the marginal-quorum reduction's r=o(1). Covariance
perturbations in Q can be amplified similarly. Neither requirement should
be silently omitted when translating the reference to fitted high-dimensional
scores. Sharper structure-aware bounds may improve these sufficient rates.

## Literature boundary

Using a heterogeneity Q statistic and noncentral chi-square inversion is
established methodology; see Dominguez Islas and Rice (2018),
[Addressing the estimation of standard errors in fixed effects meta-analysis](https://onlinelibrary.wiley.com/doi/10.1002/sim.7625),
Sections2.2 and4.1. That paper studies a weighted average of heterogeneous
effects and heterogeneity inference. The validity-count bias inequality,
its connection to a shared target estimand, arm-specific borrowing and
growing-K adaptation are the parts requiring our own full comparison and
proof. This limited search does not establish originality.

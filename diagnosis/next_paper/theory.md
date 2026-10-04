# Weak-source departures and growing K: research reference constructions

Status: 2026-09-21. These are research notes and a checked reference
construction, not a claimed new paper contribution or an established
extension of the full RoCE estimator. The current paper's frozen experiments
continue separately. The user has authorized development of this next-paper
direction.

## Target and observation model

The target estimand is a fixed, identified scalar theta, initially TATE.
There is a trusted target estimator Z0 and K source-assisted estimators Zj.
In the Gaussian reference experiment,

    (Z0, Z1, ..., ZK) = theta * 1 + (0, b1, ..., bK) + Gaussian noise.

The covariance is known in this experiment. It can include strong dependence
caused by shared target observations. At least q of the source candidates
have bj=0; their identities are unknown. Other biases are unrestricted and
may be local, large, positive or negative. K, q and the confidence-set rule
are fixed before evaluating these summaries.

The target anchor identifies theta without requiring any source to be valid.
The lower bound q is an additional borrowing assumption, not an assumption
needed for target-only identification. It must not be estimated by counting
apparently compatible sources on the same inference data and then treated as
known. Whether q can be weakened or adaptively learned requires further work.

There is a useful restriction before designing an estimator. In the
unrestricted Gaussian model with arbitrary b, write the mean as T beta,
where beta=(theta,b1,...,bK) and T has first column1 and remaining columns
(0,I_K). The first row of T^{-1} is (1,0,...,0), so the theta entry of the
inverse Fisher information is Sigma[0,0]. Thus the regular information bound
in that unrestricted model is the target-anchor variance. For linear
estimators this is elementary: unbiasedness of sum_j wj Zj for every b
requires every source weight to be zero and the anchor weight to be one.
This does not prove that every nonregular confidence procedure must equal
the target interval, but it prevents promising unrestricted first-order
variance improvement without additional structure.

Similarly, a source-weighted linear expansion with deterministic nonzero
source weights acquires a local mean shift when an otherwise unknown source
has bias h_j/sqrt(n). A confidence-set method can account for that ambiguity
without claiming a centered normal limit for its interval midpoint. This is
why coverage and interval length are the primary objectives of this prototype.

## Reference confidence set with arbitrary dependence

Choose r in {1,...,q} and allocate alpha=alpha_A+alpha_S. Let I0 be an anchor
interval with failure probability at most alpha_A. Construct source intervals
Ij such that, for every valid source,

    P(theta not in Ij) <= p,   p = alpha_S * (q-r+1) / q.

For exactly normal standardized marginal errors, Ij uses the two-sided
normal critical value Phi^{-1}(1-p/2). Define

    S_r = {t: at least r of the source intervals contain t},
    C = I0 intersect S_r.

S_r is computed as a union of closed intervals by a sweep over their
endpoints; no numerical grid or subset enumeration is required. Tied
endpoints and isolated accepted points are retained. Report the convex hull
of C if an ordinary interval is desired, retaining the original set as a
diagnostic. If C is empty, the prototype reports I0 and records the fallback.
Hulling or replacing an empty set cannot reduce coverage of the original C.

**Proposition 1 (reference coverage).** If at least q candidates are valid
and the stated marginal coverage bounds hold, then P(theta in C) >= 1-alpha.
No independence assumption, separation condition or correct classification
of the biased candidates is needed.

**Proof.** Fix any q valid candidates Vq; this is a population set, not a
selected set. Let F=sum_{j in Vq} 1{theta not in Ij}. If theta is outside S_r,
fewer than r of even these q valid intervals cover theta, so F>=q-r+1.
Markov's inequality gives

    P(theta not in S_r) <= E(F)/(q-r+1) <= q*p/(q-r+1) = alpha_S.

A union bound with the anchor failure event establishes the claim. Invalid
sources can cast arbitrary additional votes and cannot invalidate this
argument. The guarantee is uniform over their biases whenever the marginal
bounds for the anchor and valid candidates hold uniformly. QED.

The construction is a deliberately conservative benchmark. It is not claimed
to be novel, optimal, or shorter than target-only inference. A credible
comparison with Guo's searching/sampling method and RIFL is required before
making contribution claims.

## A direct growing-K implication

Suppose the anchor failure bound has an approximation error epsilon_0,n and
each valid source's failure bound has an error epsilon_j,n. The same proof
gives

    P(theta in C_n) >= 1-alpha-epsilon_0,n
                         - sum_{j in Vq} epsilon_j,n/(q_n-r_n+1).

Consequently a sufficient condition is

    epsilon_0,n -> 0,
    q_n * max_j epsilon_j,n / (q_n-r_n+1) -> 0.

If q_n/K_n -> kappa and r_n/K_n -> nu with 0<nu<kappa, the ratio q_n/(q_n-r_n+1)
remains bounded. This avoids an explicit K_n multiplier in this conservative
coverage reduction. It does NOT remove K_n from the work needed to establish
the marginal approximation errors for estimated high-dimensional nuisances.
If r_n=q_n, the bound instead requires q_n*max_j epsilon_j,n -> 0.

For an unknown majority of compatible sources, choosing
1-kappa < nu < kappa also prevents the invalid sources alone from creating
an accepted component. This is useful for studying length, but is not
required for the above coverage statement. Without further conditions there
is no asserted efficiency or oracle-length guarantee.

The current fixed-K theorem's positive limiting site fractions cannot be
maintained uniformly over infinitely many sites. A full extension needs
explicit target and source sample sizes, n_min, and their relation to K_n,
sparsity and p. The prototype does not establish those nuisance-rate results.

## Known Gaussian dependence can sharpen the reference calibration

Suppose every q-dimensional subvector of standardized valid-source errors
has the known equicorrelation matrix with correlation omega in [0,1]. Write

    Ej = sqrt(omega) U + sqrt(1-omega) Vj,

where U and the Vj are independent standard normals. Instead of the marginal
Markov bound, choose c so that

    P(the r-th smallest of |E1|,...,|Eq| <= c) = 1-alpha_S.

The same confidence-set definition then has coverage at least 1-alpha,
because these q valid candidates alone cast at least r votes with the stated
probability. Dependence between the target anchor and these candidates is
still handled by the alpha allocation and union bound.

For omega=0, if u is the (1-alpha_S) quantile of Beta(r,q-r+1),

    c = Phi^{-1}((1+u)/2).

For omega in (0,1), the CDF is the one-dimensional integral

    integral P{Binomial(q, pi(c,u)) >= r} phi(u) du,
    pi(c,u) = Phi((c-sqrt(omega)*u)/sqrt(1-omega))
               - Phi((-c-sqrt(omega)*u)/sqrt(1-omega)).

The implementation uses deterministic integration and root solving, checked
against the exact independent and fully correlated limits. For omega=1,
c=Phi^{-1}(1-alpha_S/2). This construction is currently restricted to known
Gaussian correlation. Plugging in an estimated arbitrary covariance matrix
does not inherit the guarantee without additional analysis.

For independent errors and r_n/q_n -> gamma in (0,1), the critical value
converges to Phi^{-1}((1+gamma)/2). This observation describes the calibration;
it does not by itself give an optimal confidence-interval length for the
contaminated-source model.

## Validation design and explicit limits

The Gaussian study uses K=4,8,16,32,64, per-site n=1000 and known theta=.25.
Shared standardized correlations are 0,.25,.5,.9. Invalid-source biases are
delta times the source-minus-anchor standard deviation, with delta=0,.5,1,2,4;
additional cells alternate bias signs. Thus the weak biases would scale as
n^{-1/2} along a sequence with fixed delta. These are not the current paper's
fixed-rho bounded-DGP experiments.

The assumed valid minimum is ceil(.75*K), and the vote threshold is
floor(K/2)+1. Deliberate violations with only25% unshifted sources are labeled
outside this guarantee. Comparators are target-only, pooled GLS, oracle GLS,
an illustrative naive screening-plus-GLS interval, the marginal reference,
known-correlation calibration, and deliberately incorrect independence
calibration. The naive procedure is NOT labeled as the production RoCE rule.

Known variances and Gaussian errors remove nuisance fitting from this first
experiment. Interval-center RMSE is descriptive, not a proved estimator
optimality statement. Covariance-aware studentization with estimated
covariance, arm-specific borrowing, and real three-level fitted scores remain
required parts of the full research objective.

The follow-up `validity_fraction` profile evaluates guaranteed valid fractions
of 25%, 50% and 75%, with partial versus all-valid quorum rules. It has 232
procedure settings on 100 paired data-generating settings, each with 1,000
repeats. Thus the original 75% convention is not a requirement of Proposition 1.
The variance-region prototype in variance_uncertainty.md supplies a sufficient
Gaussian estimated-variance construction; fitted-score variance regions remain
open. The two lower-bound notes distinguish coverage from uniform oracle
adaptation, and uniform_calibrated_score_remainder.md states a sufficient
conditional route to a full high-dimensional expansion.

## Literature and next proof tasks

- Guo, Kang, Cai and Small (2018), JRSSB: TSHT with voting is a selection-based
  reference. Author paper: https://zijguo.github.io/arxiv/1603.05224v3.pdf
- Guo (2023), JRSSB: searching/sampling inference addresses validity-selection
  errors under its finite-sample majority/plurality conditions. This motivates
  test inversion here, but its IV assumptions are not automatically source
  aggregation assumptions. https://doi.org/10.1093/jrsssb/qkad049
- Guo, Li, Han and Cai (2025), JASA: RIFL already treats selection uncertainty
  in federated meta-learning. Its target is a prevailing model. Our designated
  target anchor, shared-target score dependence and proposed growing-K regime
  require an explicit comparison before any novelty claim.
  https://doi.org/10.1080/01621459.2024.2443246

Outstanding: general estimated covariance; less conservative searching or
sampling calibration; validity-count sensitivity; uniform score expansions
with growing K and p; arm-specific candidate validity and joint TATE
projection; interval-length and efficiency limits; implementation using only
the allowed federated messages. The complete extension is not yet proved.

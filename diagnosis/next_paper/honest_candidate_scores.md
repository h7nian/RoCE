# Fitted candidates with a common target prediction

This is a research interface, outside the RoCE package and all current-paper
campaigns. Its purpose is to connect the shared-target score calculations to
actual fitted source models. It does not yet provide a validated confidence
interval or a replacement estimator.

## Construction and exact identities

Each site has an independent training/evaluation split, selected from its row
count and a fixed seed without treatment or outcome stratification. Given all
training observations, evaluation observations remain independent within and
between sites under the iid DGP. Source and target nuisances use the existing
fold-summed calibration kernels inside the training sample. The auxiliary
target OR minimizes an inverse-propensity weighted logistic loss; it uses
1/pi, not the odds weights of the Hou target-anchor OR.

For one arm, write m_j for the fitted source OR, q_j for its merged weight,
m_0 for this auxiliary target projection and P_tr/P_ev for target training/
evaluation averages. Define

    d_j = P_tr(m_j - m_0),
    R_ji = I(A_ji=a) q_j(X_ji) {Y_ji-m_j(X_ji)},
    mu_j_original = P_ev m_j + mean_source_eval R_j,
    mu_j_common   = P_ev m_0 + mean_source_eval (R_j+d_j).

Exactly, on every dataset,

    mu_j_common - mu_j_original = (P_tr-P_ev)(m_j-m_0).

All source candidates now share the same target evaluation prediction, while
the trusted target anchor retains its own calibrated AIPW score. The two arms
have separate common predictions. Conditional on training, source residuals
from different sites are independent; the two residual arms within a source
are generally correlated. The interface retains their full 2-by-2 covariance.
It must not be passed to a Gaussian implementation that assumes zero private
cross-arm covariance without checking that assumption.

Source summaries transform exactly. If S=sum R, T=sum R^2, and n is the
source evaluation count, the shifted summaries are S+nd and T+2dS+nd^2.
The two-arm cross sum becomes

    C + d_1 S_0 + d_0 S_1 + n d_1 d_0.

Only the source OR coefficients and existing residual moment summaries are
needed at the target. Computing d_j and the common prediction is target-local
after the model messages arrive, so this algebra does not require another
target/source exchange. The individual-level implementation is a diagnostic
reference for those summary identities.

## Why a new target projection is needed

In the correct-joint-weight / misspecified-OR branch, the source OR population
loss projects onto the target covariate distribution. Ordinary target arm OR
instead projects onto the treatment-selected target distribution; the Hou
anchor OR uses another weighting. Their pseudo-true coefficients need not
agree. The independent population audit already found a nonzero difference.

If the target propensity is correct, target inverse-propensity OR has the same
population projection as a valid source OR with correct merged weights,
assuming the same features, loss, inactive truncation at the relevant limit,
and a unique minimizer. In the correct-OR branch, both converge to the true
conditional mean under suitable nuisance consistency even if weights are
misspecified. These are assumptions about the nuisance program; the interface
does not prove them. Initial and final OR limits need not all be identical.

## Sufficient remainder conditions, not a completed growing-K theorem

Let g_j=m_j-m_0. For a valid source with a common OR probability limit,

    (P_tr-P_ev) g_j = (P_tr-P) g_j - (P_ev-P) g_j.

The evaluation term is conditionally centered and has variance bounded by
P(g_j^2)/n_ev. The training term is not conditionally centered, since its
functions were fitted on those same observations. An independent-evaluation
argument cannot bound that term.

A sufficient fixed-K condition for first-order equivalence is

    max_valid_j |(P_tr-P)g_j| = o_p(n^-1/2),
    max_valid_j ||g_j||_(L2(P)) = o_p(1),

with comparable training/evaluation sizes and suitable tails. A sparse GLM
route can bound the training term through empirical gradients and quadratic
prediction errors at order s log(pK)/n, provided the corresponding uniform
L1/L2 and empirical-process bounds actually hold. CV tuning, active truncation,
and misspecification still need explicit treatment in that argument.

For growing K, fixed-K equivalence is insufficient. To preserve an interval
scale a_nK, all original score remainders, these replacement remainders and
covariance calibration errors must be negligible relative to a_nK uniformly
over the promised valid-source class. At the illustrative dispersion scale
a_nK=n^-1/2 K^-1/4, the stated training bound would require
s log(pK) K^1/4/sqrt(n) -> 0. Additional uniform control of the evaluation
term is also needed. This condition alone neither proves that interval rate
for RoCE nor removes a nonvanishing shared-target variance floor.

For invalid sources the two OR limits can differ. The replacement can then
change their first-order behavior. They must be treated as arbitrary biased
candidates conditional on training, not as asymptotically equivalent versions
of the original invalid-source estimator.

## Remaining inference and numerical checks

The reported covariance is the empirical centered-score convention
crossprod(centered)/n^2. It is not an upper confidence bound for covariance;
even conditional iid sampling would require n/(n-1) to make this particular
sample covariance estimator unbiased. Finite-sample variance calibration and
valid-candidate conditional bias bounds remain open.

The numerical radius12 only gives a very large worst-case residual bound
through exp(12). Substituting that bound into a bounded-score reference can
produce uninformative intervals. Known covariate support and the actual fitted
coefficients may give tighter conditional bounds, but observed sample maxima
are not valid population bounds by themselves.

The first implementation check used C2, n=1000/site, p=4, one source,
three calibration blocks and100 nuisance lambdas. All29 assertions passed:
evaluation mutation preserves the split, hashes, target/source/common fits
and training shifts; training mutation is rejected; independent IPW KKT
residuals are below1e-11; summary and covariance reconstructions agree.
Artifacts: scratch implementation/next_paper/v9/honest_candidate_checks_v1/.
C3 with two sources/two-round passed37 assertions, including evaluation
mutation followed by fitting in an independent new cache. C1 with two sources/
one-round passed41 assertions, also comparing proximal Newton and coordinate
descent for the auxiliary loss at the same penalty. The coefficient comparison
used tolerance1e-6 and the objective difference was below1e-9. These checks are
under honest_candidate_checks_v3/ and v4/; the first two-source test in v2/
failed because its fixture passed a scalar source sample-size vector, which
explicitly requests one source. That failure is retained. The corrected runner
uses a length-K vector and checks the realized count before fitting. Production
method/baseline workers already pass that vector; sampled K4/6 production
outputs have total sizes5000/7000 as intended.

The later [independent covariance pilot](honest_covariance_pilot.md) adds an
optional held-out index list without changing the default packet. It separates
covariance estimation from final means and checks both the site boundaries and
unequal-sample-size scaling. This is research-only and leaves production
cross-fitting unchanged.

Broader dimensions, validity patterns and ultimately uniform bias/coverage
validation are still required. Current-paper jobs remain unchanged.

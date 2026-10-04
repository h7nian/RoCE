# Valid fitted-score drift has two different components

Condition on an honest training sample. Let z contain the target anchor and
source candidates for both arms, with conditional mean H theta + r and
conditional covariance Sigma. For a truly valid coordinate subset V, let
H_V repeat the two arm means and assume Sigma_V is positive definite.
Validity here refers to the data-generating transportability condition; it
does not assert that fitted-score drift r_V is exactly zero.

Write P = inverse(Sigma_V), I = H_V' P H_V and define

    delta = inverse(I) H_V' P r_V,
    r_perp = r_V - H_V delta.

Then H_V' P r_perp = 0 and the exact Pythagorean identity is

    r_V' P r_V = delta' I delta + r_perp' P r_perp.

The oracle-valid GLS TATE estimator has conditional mean drift
delta_1 - delta_0. In contrast, its Gaussian residual dispersion statistic
only has noncentrality Lambda = r_perp' P r_perp, with |V|-2 degrees of
freedom. Adding the same arm-specific shift H_V a to every valid candidate
changes delta to delta+a but leaves Lambda unchanged.

Thus small candidate disagreement alone cannot establish small fitted-score
bias relative to the causal estimand. A common location shift is invisible
to residual dispersion. This is an algebraic limitation of these summaries;
it is not a claim that the original causal observed-data model is unidentified.
The fitted method still needs a bound or negligible-rate argument for the
valid anchors and candidate nuisance drift. Covariance checks cannot replace it.

## A growing-K sensitivity calculation

Suppose |V| is proportional to K, the smallest eigenvalue of Sigma_V is at
least c/n, and max |r_V| is bounded by b_n. Then

    Lambda <= r_V' inverse(Sigma_V) r_V <= C n K b_n^2.

For this sufficient upper bound to be negligible relative to the central
Gaussian residual fluctuation of order sqrt(K), it suffices that

    sqrt(n) K^(1/4) b_n -> 0.

This is a sufficient requirement for ignoring valid-drift contamination of
that dispersion statistic at its central fluctuation scale. It is not a
necessary condition for overall interval coverage, and a noncentrality-based
procedure could explicitly allow contamination instead. A root-n nuisance
remainder alone does not imply this stronger condition as K grows. If a proved
b_n is of order s log(pK)/n, the corresponding sufficient restriction is
s log(pK) K^(1/4) = o(sqrt(n)). CV-selected fits do not automatically satisfy it.

This requirement also says nothing by itself about the common mean component
delta: that component needs separate control relative to the final estimator's
standard deviation. Different covariance eigenvalue regimes require their own
scaling; the shared target component must be retained in that calculation.

For example, in a scalar shared-component model with oracle variance
v_t/n + v_s/(nK), negligible TATE drift d_n requires

    sqrt(n) |d_n| / sqrt(v_t + v_s/K) -> 0.

If v_t is bounded below, ordinary root-n drift control suffices for this
particular requirement. If v_t=0 and v_s is bounded away from zero and infinity,
the condition becomes sqrt(nK) |d_n| -> 0. A second-order nuisance drift of
order s log(pK)/n would then require s log(pK) sqrt(K) = o(sqrt(n)). The latter
is stronger than the quarter-power condition for residual dispersion. This is
an explicit covariance-model example, not a universal RoCE variance formula.
It explains why the fitted diagnostic reports drift relative to its actual
oracle SD as well as residual dispersion: adding sources can shrink the
relevant SD while leaving a shared nuisance error insufficiently controlled.

## Diagnostic use

honest_bias_decomposition.R computes both components from actual fitted
conditional-population packets and oracle validity labels. Tests add a common
arm shift and check invariance of residual dispersion and the exact change in
TATE drift. No Gaussian samples or proposed confidence intervals are needed
for this algebra check. In the bounded observed-data study, the Mahalanobis
quantity is a diagnostic; it is not declared exactly chi-square distributed.

This complements the conditional mean/covariance study. It does not supply a
feasible nuisance-bias envelope, a fitted covariance uncertainty region or a
completed growing-K theorem.

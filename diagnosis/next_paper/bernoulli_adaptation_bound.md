# A direct observed-data local-adaptation lower bound

This transfers the earlier Gaussian lower-bound argument to a simple causal
observed-data submodel. It does not establish an upper bound for the full
high-dimensional method or claim originality.

At every site observe n independent records. Treatment A is randomized with
known probability e in (0,1). Conditional on A=1, Y is Bernoulli with site
mean m_j; the control outcome law is identical at all sites and unchanged
throughout the argument. Covariates may have any common bounded distribution
independent of these variables. Potential outcomes can be generated
independently of treatment, so consistency, overlap and unconfoundedness hold.
The target treatment effect changes exactly with its treated mean. At least
q=ceil(gamma K) source treated means agree with the target; control sources
are all valid. Constant logistic outcome and weight models put this submodel
inside a sparse high-dimensional model class.

## Exact likelihood identity

Let the all-valid null have treated mean m at every site. For one observed
record, the likelihood ratio L_h for changing that mean to m+h equals one
when A=0 and is the Bernoulli likelihood ratio when A=1. Therefore

    E_m[L_h L_g] = 1 + e h g / {m(1-m)}.

For n independent records, the cross-moment is the nth power. Define
sigma_n^2=m(1-m)/(ne), u=delta/sigma_n, and choose w in (gamma,1), with
c=w/(1-w). Under an alternative, the target mean is m+delta; each source mean
is m+delta with probability w and m-c delta otherwise. Both probabilities
remain in (0,1) for sufficiently small delta. The mixture's source mean equals
m, but its valid sources agree with the alternative target m+delta.

For one source the exact chi-square divergence from the null is

    D_n(u) = w^2(1+u^2/n)^n
             + 2w(1-w)(1-c u^2/n)^n
             + (1-w)^2(1+c^2 u^2/n)^n - 1.

There is a useful nonnegative expansion for stable numerical evaluation:

    D_n(u) = sum_(r=2)^n choose(n,r) (u^2/n)^r
                         [w+(1-w)(-c)^r]^2.

The r=1 term vanishes. For n>=2 and small u,

    D_n(u) = (1-1/n)c^2 u^4/2 + O(u^6).

The source mixture for n=1 is exactly the null Bernoulli law; no Gaussian
approximation is being invoked in any of these identities.

## Valid-count conditioning and length consequence

The number of valid sources under this prior is Binomial(K,w). Condition on
it being at least q. Every resulting parameter vector satisfies the declared
valid-count condition. The conditioning cost in total variation is at most
epsilon_K=P{Binomial(K,w)<q}, which decays exponentially for fixed w>gamma.

Before conditioning, site independence and the target likelihood ratio give

    chi2(Q,P0) = (1+u^2/n)^n [1+D_n(u)]^K - 1,
    TV(P0,Q_valid) <= 0.5 sqrt(chi2(Q,P0)) + epsilon_K.

Set u=a K^(-1/4), with sufficiently small fixed a>0. For n>=2, including
joint sequences with n increasing, the divergence remains bounded by a small
constant. A confidence interval uniformly covering the target effect with
probability at least 1-alpha must contain both the null and alternative
effects under P0 with probability at least 1-2alpha-TV(P0,Q_valid). Hence

    E_0[length(I)] >= c_(gamma,alpha,m,e) / {sqrt(n) K^(1/4)}.

Here the oracle knows which sources are valid, not the unknown outcome
probabilities themselves. This is an observed-data local-adaptation obstruction, not population
nonidentifiability. The trusted target identifies its effect. It rules out
uniformly attaining the all-valid oracle order 1/sqrt(nK) while retaining
coverage over the unrestricted local-bias class. It does not forbid an oracle
result in a separate, strongly separated regime.

The exact likelihood calculations were checked by enumerating all four
observed (A,Y) categories for n=1,...,4; ten checks passed. Finite positive
lower bounds for n=1000, K=4,...,1024 and majority fractions .6/.75 are saved
in scratch v9/bernoulli_adaptation_bound/finite_bounds.csv.

The dispersion expected-length upper bound in dispersion_expected_length.md
currently holds in the known-Gaussian reference. Matching that order with
observed data, estimated covariance and high-dimensional fitted nuisances
remains a separate task; the lower bound alone does not validate our estimator.

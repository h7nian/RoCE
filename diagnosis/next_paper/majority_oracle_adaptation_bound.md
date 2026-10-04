# Local oracle adaptation can fail even with a guaranteed majority

Status: a lower bound in the independent Gaussian-summary model. It does not
establish a minimax-optimal rate or an observed-data RoCE lower bound. No
novelty claim is made; its purpose is to set a defensible research objective.

See `adaptive_gaussian_literature_update.md` before making any novelty claim.
An exact binary observed-data extension is in `bernoulli_adaptation_bound.md`;
it does not establish an upper bound for the full fitted method.

## Parameter class

Observe independent Z_0,...,Z_K with common known variance sigma_n^2=1/n.
The target anchor has mean theta. At least q_K=ceil(gamma*K) source means
equal theta, where 0<gamma<1 is fixed. Other source means are unrestricted.
Let I be an ordinary interval with coverage at least 1-alpha uniformly over
this class, with alpha<1/2.

Even at the all-valid law P_0 with theta=0 and all source means zero, there is
a positive constant c(gamma,alpha) such that, for sufficiently large K,

    E_0[length(I)] >= c(gamma,alpha) / (sqrt(n)*K^(1/4)).

The all-valid oracle interval has length of order 1/sqrt(n*K). Thus uniform
coverage over the unrestricted local-bias class precludes uniformly attaining
that oracle length, including adaptation at the all-valid law. This does not
preclude oracle behavior under additional separation or structural conditions.

## A prior supported on local alternatives after conditioning

Choose pi in (gamma,1), for example pi=(1+gamma)/2, and put c=pi/(1-pi).
Let the alternative target be theta=delta>0. Independently for each source,
draw its mean to be delta with probability pi and -c*delta otherwise. The
number of valid sources is B~Binomial(K,pi). Condition this prior on B>=q_K.
Every parameter vector in the conditioned prior belongs to the assumed class.
Write Q for the unconditioned mixture law and Q_valid for the conditioned law.

The conditioning cost in total variation is bounded by

    TV(Q,Q_valid) <= P(B<q_K) =: epsilon_K,

which tends to zero exponentially because pi>gamma. All alternative parameter
vectors share the same target theta=delta, so uniform coverage also holds
after averaging over the conditioned prior.

## Exact chi-square calculation

Set u=delta/sigma_n. For one source, the standardized mixture has density
pi*phi(x-u)+(1-pi)*phi(x+c*u), whereas the null density is phi(x). Its mean
is zero. The chi-square divergence from the null is exactly

    D(u) = pi^2*exp(u^2) + 2*pi*(1-pi)*exp(-c*u^2)
           + (1-pi)^2*exp(c^2*u^2) - 1
         = c^2*u^4/2 + O(u^6).

The quadratic term cancels because the mixture mean is zero. Independence
and the anchor mean shift give

    chi2(Q,P_0) = exp(u^2) * (1+D(u))^K - 1,
    TV(P_0,Q_valid) <= .5*sqrt(chi2(Q,P_0)) + epsilon_K.

Choose u=a*K^(-1/4). The chi-square divergence converges to
exp(c^2*a^4/2)-1. Taking a>0 sufficiently small makes the resulting TV upper
bound strictly less than 1-2*alpha. For example a=A/sqrt(c) yields limiting
bound .5*sqrt(exp(A^4/2)-1); A can be chosen independently of gamma.

## Interval-length conclusion

Coverage at P_0 and under Q_valid implies

    P_0({0,delta} subset I) >= 1-2*alpha-TV(P_0,Q_valid).

Multiplying this positive lower bound by delta=a*sigma_n*K^(-1/4) proves
the expected-length result. Conditioning the prior is essential: a Bernoulli
mixture without that step would occasionally violate the guaranteed valid count.
For disconnected confidence sets the conclusion concerns diameter/hull rather
than Lebesgue measure.

## Interpretation and limits

This is a uniform local-efficiency obstruction, not population
nonidentifiability. The target is always identified by its trusted anchor.
Alternative bad-source biases are -(c+1)*delta, smaller than the per-site
standard error as K grows. These alternatives are not excluded by a fixed
positive valid fraction, even one exceeding one half.

The K^(-1/4) lower bound is not claimed sharp for the whole class; other least
favorable priors can potentially be harder. It also does not say our current
vote confidence sets attain this rate. The known-covariance experiments must
report their actual length, not describe nominal coverage as oracle efficiency.
When gamma<=1/2, the separate two-point bound in half_valid_lower_bound.md is
stronger; this mixture construction is useful because it also covers gamma>1/2.

A shared target component can give the oracle a variance floor, so comparisons
with 1/sqrt(n*K) apply to the independent submodel here. Transferring the result
to calibrated source-assisted TATEs requires embedding the local alternatives
in admissible observed-data models and controlling the approximation error.
RIFL's separated-regime oracle theorem and this local-bias bound address
different parameter regimes.

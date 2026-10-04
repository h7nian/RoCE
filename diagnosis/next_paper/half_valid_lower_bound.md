# A local ambiguity bound when at most half the sources are guaranteed valid

Status: a Gaussian-summary lower bound derived for the proposed research
model. Embedding these alternatives in the full calibrated RoCE experiment
remains to be proved. This is not a population nonidentifiability claim.

## Model and claim

Observe Z ~ N(mu, Sigma), with known positive-definite Sigma of size K+1.
Index 0 is a trusted target anchor: mu_0=theta. At least q source means equal
theta; all other means are unrestricted. Suppose 1 <= q <= floor(K/2).
Write sigma_0^2=Sigma_00. A reported ordinary interval I must cover theta
with probability at least 1-alpha at every distribution in this class.

For any t>0 there are two distributions P_0,P_1 in the class with target
parameters 0 and delta=t*sigma_0 such that

    TV(P_0,P_1) = 2*Phi(t/2)-1,
    E_0[length(I)] >= t*sigma_0 * [2-2*alpha-2*Phi(t/2)]_+.

For alpha=.05 and t=1 the lower bound is approximately .5171*sigma_0.
The same statement holds under P_1. It is an order bound, not a sharp
constant or a proof that the ordinary target-only interval is optimal.
For a disconnected confidence set the argument bounds its diameter/hull,
not its Lebesgue measure.

## Construction and proof

Choose disjoint source index sets A,B of size q. Let

    c_j = Sigma_j0 / Sigma_00,
    Delta = delta * Sigma e_0 / Sigma_00.

Define mu^(0)_0=0 and mu^(1)_0=delta. For source coordinates set

    j in A: mu^(0)_j=0,             mu^(1)_j=c_j*delta;
    j in B: mu^(0)_j=(1-c_j)*delta, mu^(1)_j=delta;
    others: mu^(0)_j=0,             mu^(1)_j=c_j*delta.

Each model has at least q valid sources, although their guaranteed-valid
sets differ. Their mean difference is exactly Delta. Hence

    Delta' Sigma^{-1} Delta = delta^2 / Sigma_00 = t^2.

Equal-covariance normal laws at Mahalanobis distance t have total variation
2*Phi(t/2)-1. This follows directly by integrating the likelihood-ratio
test over the half-space closer to mu^(1); its two error probabilities
are both Phi(-t/2). Also KL(P_0,P_1)=t^2/2.

Coverage under P_1 gives P_0(delta in I) >= 1-alpha-TV(P_0,P_1). Combining
this with P_0(0 in I) >= 1-alpha yields

    P_0({0,delta} subset I) >= 1-2*alpha-TV(P_0,P_1).

An interval containing both points has length at least delta, proving the
claim. Exchanging the two models proves the symmetric statement.

## Growing-K interpretation

For independent equal-variance summaries, Sigma=I/n, the construction has
zero source means on A and delta on B under both models. The source data
alone cannot distinguish these two mean labels, while the target anchor
does distinguish them with fixed nonzero testing error at delta=t/sqrt(n).
Thus an interval uniformly valid over these local alternatives must have
expected length of order at least 1/sqrt(n), irrespective of K.

If the q valid sources were known, their average with the target would have
standard error 1/sqrt(n*(q+1)). Consequently, when q increases with K but
q <= K/2, uniformly attaining that oracle interval-length rate is impossible
in this Gaussian class without additional restrictions. Pointwise gains
under stronger separation remain possible.

With a shared target component, the covariance-aware construction above is
essential: keeping all source means identical across the two alternatives
would instead give a conditional-anchor information bound and could miss
the least favorable direction. The general bound uses the marginal anchor
variance Sigma_00. An oracle also faces the common-noise variance floor, so
its width need not shrink like 1/sqrt(n*K) in this correlated case.

The anchor still identifies theta in the population under both models.
The result concerns uniform finite-sample/local efficiency with unknown
source labels, not failure of population identification. It does not extend
unchanged to q>K/2, because disjoint valid sets of size q cannot then exist.

## Consequences for the research design

- Evaluate coverage separately from oracle-relative length. A method can
  remain valid without matching oracle efficiency near this boundary.
- Include 25%, 50% and 75% guaranteed-valid fractions, with a valid quorum
  r<=q specified in advance. The reference coverage proof does not require
  q>K/2; preventing an invalid group from forming its own accepted component
  is a separate length/identification-of-source-cluster issue.
- The current-paper half-invalid K=4/6 controls are useful boundary examples.
  Their fixed rho is not the local-alternative sequence used in this proof.
- To claim a full RoCE impossibility result, construct corresponding observed
  data laws with the required transport/calibration assumptions and justify
  the Gaussian approximation uniformly. This note does not supply that step.

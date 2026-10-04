# From independent site scores to the dispersion reference

This is a conditional approximation lemma for the honest fitted-score design.
It is not a completed theorem for the existing full cross-fitted estimator.
The derivation below combines a published fixed-dimensional CLT with the
shared-target score structure; no new general CLT rate is claimed.

## Conditions and block representation

Condition on the training data and all fitting randomness. Evaluation rows
must be independently sampled from training and from other sites. The honest
500/500 split uses index-only randomization for this reason. Outcome-stratified
splits or conditioning on every cross-fitting training fold do not automatically
supply this conditional independence.

The target supplies a block of at most four scores: its two anchors and two
common predictions. Each source supplies two residual scores. Center and
standardize each block on its covariance range. Let

    U_j = n_j^(-1/2) sum_i V_ji,
    E[V_ji]=0, Cov(V_ji)=I_(d_j), d_target<=4, d_source<=2.

Blocks are independent, and observations are iid within a block. Require a
uniform conditional fourth-moment bound E||V_ji||^4<=M. The nuisance feature
dimension p can grow; it enters the fitting/bias and moment conditions, not
the block dimensions d_j. Bounded covariates alone do not prove a uniform
bound after score standardization.

The complete candidate vector has the exact conditional representation

    Z = H theta + r + L U,    L L' = Sigma.

Training shifts are part of r. Private cross-arm covariance is retained in L.
Assume Sigma is positive definite on the candidate space used by the procedure.

## Published coupling input

For centered iid vectors with identity covariance and finite fourth moment,
Bonis (2020), Theorem1, equation8, gives

    W2(U_j,G_j) <= C d_j^(1/4)
      ||E[V_ji V_ji' ||V_ji||^2]||_HS^(1/2) / sqrt(n_j),    C<14.

The matrix norm is at most E||V_ji||^4 before taking its square root. Thus
fixed block dimensions and bounded standardized fourth moments give squared
coupling errors e_j^2=O(1/n_j). The original theorem and its assumptions were
checked on PDF page5. [Bonis, primary manuscript](https://arxiv.org/pdf/1905.13615)

Choose the couplings independently across sites. Their centered errors have
total squared norm bounded in expectation by E_n=sum_j e_j^2.

## A bound uniform over invalid-source means

For the reduced bias map M0, write V=M0 Sigma M0' and
A=V^(-1/2) M0 L. Then A A'=I_nu, so its operator norm is one. With
delta=V^(-1/2) M0 r,

    sqrt(Q)=||delta+A U||,
    R_G=||delta+A G||,
    E|sqrt(Q)-R_G|^2 <= E_n.

R_G is noncentral chi with nu>=1 degrees of freedom. Its density is bounded
by a universal constant, uniformly over nu and delta. One way to see this is
the Poisson mixture of central chi distributions. Central chi with df1 has
density at most sqrt(2/pi). For df>=2, evaluating the density at its mode and
using the ordinary Stirling lower bound for Gamma bounds it by sqrt(e/pi)<1.
The mixture therefore also has density at most1. The mixture identity is
documented in the [R statistics reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Chisquare.html).

For any t>0, coupling and Markov's inequality consequently give

    sup_x |P(Q<=x)-P(chi-square_nu(||delta||^2)<=x)|
      <= t + E_n/t^2
      <= 3*2^(-2/3) E_n^(1/3),

where the last line chooses t=(2 E_n)^(1/3), and the bound can be capped at1.
It is uniform in the mean vector r, including weak or strong invalid-source
errors; no separation condition enters this distributional step. For comparable
site sizes, E_n=O(K/n_min), giving the sufficient regime K=o(n_min).

Scalar GLS noise and reference noise admit the same coupling-to-CDF argument.
Independence between their noise and Q is unnecessary because interval coverage
uses a union bound. With the actual conditional covariance and justified
valid-score bias budgets, the Gaussian interval argument therefore transfers
to honest observed scores with an additional O((K/n_min)^(1/3)) approximation
error, plus the failure probability of the training/moment/bias conditions.
These constants need not be useful at n=1000; this is an asymptotic route.

## Mean drift remains a separate condition

The bias-budget inequalities hold for arbitrary r once the valid coordinates
satisfy their envelopes. A zero-budget implementation can only omit those
envelopes asymptotically if their induced GLS allowances are negligible relative
to the final noise SD. A sufficient direct condition is B_n/sqrt(v_all)->0,
including the reference allowance. Under bounded coefficient L1 norms this
reduces to b_n/sqrt(v_all)->0.

For a variance floor of order1/n, root-n drift control suffices. If the relevant
variance can instead be of order1/(nK), a drift b_n=O(s log(pK)/n) requires
s log(pK) sqrt(K)=o(sqrt(n)). The quarter-power condition in the earlier
decomposition note concerns ignoring drift in Q at its central fluctuation
scale; it is not by itself the condition for interval coverage. Keeping the
drift as noncentrality and retaining the separate mean allowance avoids that
confusion.

## Estimated covariance: a sufficient independent-pilot route

Suppose a covariance pilot is independent of final mean evaluation, conditional
on model training. Condition on that pilot too, and assume

    (1-epsilon) Sigma_hat <= Sigma <= (1+epsilon) Sigma_hat, epsilon<1.

Construct GLS, the reduced map and the deterministic bias bounds using
Sigma_hat. The normalized Gaussian noise covariance is then between
(1-epsilon)I and(1+epsilon)I. Coupling this noise to a standard normal adds
at most O(nu epsilon^2) to the squared error above. Thus a sufficient total
approximation error is

    O((K/n_min + K epsilon^2)^(1/3)).

If separately justified block covariance concentration yields
epsilon=O_p(sqrt(log K/n_cov)), the sufficient additional restriction is
K log K=o(n_cov). This rate needs stronger uniform boundedness or tail
conditions for standardized scores; the fourth-moment CLT alone does not
establish it. Add the covariance event's failure probability as well.

The original500/500 diagnostic estimates covariance and means from the same
evaluation rows, so it does not satisfy this independent-pilot premise. The
[independent-pilot diagnostic](honest_covariance_pilot.md) now implements and
checks500/250/250 training/pilot/evaluation, with covariance components scaled
to the final evaluation sample sizes. This checks sample boundaries and the
statistical interface, not the covariance concentration premise or coverage.
Reusing training-score covariance requires different
empirical-process arguments. No covariance-confidence construction or full
cross-fitting extension is supplied here.

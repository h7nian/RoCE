# From calibrated arm scores to an honest pairwise TATE confidence set

This is a conditional reduction, not a completed theorem for fitted RoCE.
It specifies what the nuisance expansion and variance estimation must prove.
No independence among candidate estimators or among arm pairs is required.

## Quantities and assumptions

Choose fixed, population-valid arm subsets of sizes q1+1 and q0+1, including
the trusted target in each arm. Their Cartesian product has Q=(q1+1)(q0+1)
valid contrasts. For each such pair u, write

    Dhat_u - tau = L_u + R_u,    sigma_u^2 = Var(L_u) > 0.

The anchor contrast has the analogous expansion with subscript A. Suppose
the standardized leading terms have marginal CDF errors at most epsilon_P
and epsilon_A relative to the standard normal. The CDF approximation is
unconditional; it need not hold after conditioning on fitted nuisances.

Suppose a single event E, with P(E) >= 1-beta, gives all the following:

    |R_u|/sigma_u <= r_P,
    |sigmahat_u/sigma_u - 1| <= delta_P < 1,   for every selected valid pair;
    |R_A|/sigma_A <= r_A,
    |sigmahat_A/sigma_A - 1| <= delta_A < 1.

Invalid pairs may have arbitrary biases and dependence. Their apparent
agreement is not used to declare or estimate q1 and q0.

## Coverage bound

Allocate alpha=alpha_A+alpha_S, choose r<=Q, and let

    p = alpha_S * (Q-r+1)/Q,
    c_P = Phi^{-1}(1-p/2),  c_A = Phi^{-1}(1-alpha_A/2),
    d_P = r_P+c_P*delta_P,   d_A = r_A+c_A*delta_A,
    eta_Q = Q/(Q-r+1).

The implemented intervals use Dhat_u +/- c_P*sigmahat_u. On E, failure of a
valid pair interval implies |L_u|/sigma_u > c_P-d_P. The Gaussian CDF bound
and the maximum normal density therefore give the unconditional proxy bound

    P(|L_u|/sigma_u > c_P-d_P)
      <= p + 2*epsilon_P + sqrt(2/pi)*d_P.

The same argument holds for the anchor. The inequality stays a valid upper
bound if c_P-d_P is negative, although it then need not be informative.

If fewer than r of the selected Q valid intervals cover tau, at least
Q-r+1 of their failure proxies occur. Apply Markov's inequality to that
count, then union with anchor failure and E-complement. This yields

    P(tau is in the reported interval)
      >= 1-alpha-beta
         -2*epsilon_A-sqrt(2/pi)*d_A
         -eta_Q * (2*epsilon_P+sqrt(2/pi)*d_P).

The event E-complement is charged once, not once per pair. The code's
intersection, recorded empty-set anchor fallback, and final hull cannot
reduce coverage below this bound: fallback and hull only enlarge the set
that enters the failure argument. When q1=q0=0 the procedure instead uses
the ordinary target interval, without an artificial source error allocation.

## Consequence for growing K

For partial quorum r <= nu*Q with fixed nu<1,

    eta_Q <= 1/(1-nu),
    p >= alpha_S*(1-nu).

Thus c_P and eta_Q stay bounded as K grows. The reduction itself does not
multiply the marginal approximation error by K^2. Full quorum r=Q instead
has eta_Q=Q and a growing critical value; it needs substantially stronger
uniform approximation bounds. These are distinct inference guarantees,
even if finite Gaussian experiments give similar empirical coverages.

This removes one combinatorial error multiplier. It does not prove the
uniform score expansion, Gaussian approximation or variance event E.

## What joint-arm moment estimation must control

Let the joint arm covariance be Sigma, ordered as all mu1 candidates then
all mu0 candidates. The pair variance is

    sigma_jk^2 = Sigma_(1j,1j) + Sigma_(0k,0k) - 2*Sigma_(1j,0k).

If the largest covariance-entry error is e, the largest pair-variance error
is at most 4e. Under sigma_jk^2 >= c/n, uniform relative variance consistency
therefore follows from n*e=o_p(1). The corresponding relative SE bound uses

    |sqrt(1+x)-1| = |x| / (sqrt(1+x)+1).

Cross-arm entries contain the target contribution and the same-source
residual covariance. Treating arm residuals as independent changes this
variance even when source observations are mutually independent across sites.

If each arm has a uniform absolute expansion remainder bounded by r_abs,
each pair remainder is at most 2*r_abs. With sigma_jk >= sqrt(c/n), this
gives r_P <= 2*sqrt(n/c)*r_abs. There are no new nuisance fits for K^2 pairs.

For orientation, if the earlier conditional calibration expansion establishes
r_abs=O_p(s*log(p*K*F)/n) with a uniform high-probability bound, the familiar
s*log(p*K*F)=o(sqrt(n)) condition makes the standardized remainder vanish.
This statement is conditional on that expansion; it is not a proof that the
current calibrated fits achieve the necessary uniform rates under every
misspecification, clipping boundary or growing-source regime.

The empirical covariance extractor is now numerically checked against saved
arm and TATE scores. Those identities establish correct moment assembly;
they do not establish n*e=o_p(1) or the distributional error bounds above.

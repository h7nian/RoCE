# Joint dispersion with a trusted target and arm-specific validity counts

This is a known-Gaussian reference construction and a sufficient-condition
research route. It is not a novelty claim, an estimated-covariance theorem, or
validated inference for the fitted high-dimensional RoCE estimator.

## Model and known-valid information

Stack all treated-mean candidates followed by all control-mean candidates in
X ~ N(B theta + b, Sigma). Here theta=(mu1,mu0), B contains the two arm
indicators, and the estimand is ell'theta with ell=(1,-1). Each arm contains
the trusted target followed by K source-assisted candidates. Sigma is known
and positive definite. At least q1 and q0 sources are valid for the respective
means, possibly with different identities. Biases of other candidates are
arbitrary under this Gaussian model.

Let V0 contain both targets and all sources in an arm ONLY when its declared
minimum valid count equals K. These coordinates are known valid by assumption.
Do not enlarge V0 using apparent agreement in the observed data. Let U be the
remaining coordinates.

The reference interval uses the BLUE for TATE based on V0. This exploits
exactly the declared information. In particular, when the control arm is
declared fully valid, the comparison must also pool its information rather
than use an unnecessarily noisy target-only benchmark.

## Joint mean and bias norm

Define I=B'Sigma^(-1)B and a=Sigma^(-1)B I^(-1)ell. The pooled estimate is
M=a'X, with variance v=ell'I^(-1)ell and bias a'b. Define

    P = Sigma^(-1) - Sigma^(-1)B I^(-1)B'Sigma^(-1).

For an admissible valid subset V (containing V0 and the required valid counts),
let v_V be the BLUE variance for the contrast from X_V. Then

    (a'b)^2 <= (v_V-v) b'P b, whenever b_V=0.

To see the sharp constant, put S=V-complement. P_SS is positive definite
because B_V has full column rank. Cauchy-Schwarz gives the optimal constant
a_S'P_SS^(-1)a_S. This equals v_V-v: the restricted BLUE is the least-variance
unbiased contrast supported on V, and its noise covariance with the full BLUE
is v. Thus its squared distance from the full BLUE is the variance difference.
Equivalently, this is a constrained quadratic minimization. Tests verify both
identities and equality for the bias direction P_SS^(-1)a_S.

Consequently C=max_V(v_V-v) bounds the contrast bias by sqrt(C lambda), where
lambda=b'P b. No classification of individual sources is needed.

## Remove residual noise from known-valid coordinates

The full lack-of-fit statistic has 2K residual degrees of freedom. If one arm
is fully valid by assumption, some of these residuals contain no bias signal.
Use the nested-model statistic comparing the common-mean model B theta with
the alternative B theta + E_U b_U. Its degrees of freedom are |U|, and

    Q ~ noncentral_chisquare(|U|, lambda),    lambda=b'P b.

The noncentrality is unchanged because b is supported on U and hence belongs
to the alternative model. Only noise dimensions are removed. For example,
with all control sources known valid, the degrees of freedom fall from 2K to K.

Implementation computes the GLS estimator of unrestricted b_U. If
N=Sigma_UV0 Sigma_V0V0^(-1), A is the BLUE coefficient matrix for theta from V0,
and D=B_U-N B_V0, then

    b_hat_U = X_U - N X_V0 - D A X_V0 = R X,
    Q = b_hat_U' (R Sigma R')^(-1) b_hat_U.

R B=0 and Cov(a'X,RX)=0. Its inverse covariance equals P_UU. Independent tests
also reconstruct Q as the full-model residual quadratic form minus the
known-valid residual quadratic form. If U is empty, use the ordinary full
BLUE confidence interval without a dispersion-error allocation.

## Interval and its coverage scope

Invert the noncentral chi-square distribution to get an upper confidence
bound U_delta(Q) on lambda. With noise error alpha_N and dispersion error
delta, an interval is

    M +/- [z_(1-alpha_N/2) sqrt(v) + sqrt(C U_delta(Q))].

A union bound gives coverage at least 1-alpha_N-delta uniformly over b in the
stated validity-count class. This includes arbitrarily weak nonzero biases.
If intersecting with the known-valid reference interval, allocate a separate
alpha_R and require alpha_N+delta+alpha_R<=alpha. Replacing an empty
intersection by the reference interval cannot reduce coverage of the
intersection procedure. When no source validity is guaranteed in either arm,
the implementation directly uses the target interval.

This argument requires Gaussian approximation around ALL biased candidate
means. Valid-source-only asymptotic linearity is insufficient. The earlier
one-round initialization audit remains relevant. A wrong minimum-valid-count
assumption can invalidate borrowing even though the target still identifies
TATE; the Gaussian study explicitly includes that failure case.

## Exact subset selection under shared predictions

Suppose the source covariance in each arm is H_aa 11'+diag(d_a,j), the
cross-arm source covariance is H_10 11' (zero local cross-arm residual
covariance), and each target/source cross-covariance is constant across sources
within an arm. The associated common 4-by-4 target/prediction covariance must
be positive semidefinite.

For a fixed number of valid sources in each arm, the restricted BLUE variance
depends on source identities only through their arm-specific private
precisions. Increasing any selected private variance can only increase BLUE
variance; the common covariance is unchanged. Therefore the worst admissible
subset consists of the q_a largest private variances in each arm. Because
cross-arm private covariance is zero and validity sets may differ, the two
subset choices can be made independently.

This replaces choose(K,q1)*choose(K,q0) subset enumeration with O(K log K)
subset selection. The current dense covariance checks/GLS operations still
cost O(K^3); the entire implementation is not claimed to be O(K log K).
The structural option verifies the covariance equalities and rejects matrices
outside the structure. It never silently projects an empirical covariance onto
this form. A general exact-enumeration option and a known-valid upper bound
remain explicit alternatives.

The same selection rule gives an exact fast comparator for the earlier
armwise dispersion intervals. The Gaussian comparison keeps the error budget
and known-valid reference the same when isolating joint versus armwise bias
control. The prior ordinary-target-reference version is retained separately.

## Relation to the literature and outstanding work

Guo et al. (2018) study TSHT with voting and a plurality identification rule;
Guo (2023) addresses inference errors caused by locally invalid IVs using
searching and sampling. These motivate controlling uncertainty from source
validity rather than relying on perfect selection. The target here is already
identified by its own data; the validity-count assumption governs borrowing.
Our quadratic construction is a different Gaussian reduction, not a direct
application of their IV theorems or a claim to improve their methods.

Primary sources checked: [Guo et al. (2018)](https://academic.oup.com/jrsssb/article/80/4/793/7048331),
[Guo (2023)](https://academic.oup.com/jrsssb/article/85/3/959/7174915).
The latter also discusses an oracle bias-aware comparison, which reinforces
the need to distinguish oracle selection and bias-aware efficiency benchmarks.

Remaining requirements include estimated-covariance uncertainty, uniform
nuisance expansions, growing-K quadratic-form approximation, expected length,
and actual fitted n=1000 experiments. At OR-correct limits, common source
predictions and zero private cross-arm residual covariance are natural; local
departures may preserve that limiting structure. Finite fitted predictions
are not identical, and fixed heterogeneous sources need not have this
structure. These links must be proved or explicitly bounded before claiming
the fast Gaussian reference applies to the full method.

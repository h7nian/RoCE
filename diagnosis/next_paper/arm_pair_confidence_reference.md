# Arm-specific borrowing through a pairwise confidence-set reference

Status: the reference construction and joint-arm moment interface are now
implemented and tested. The conditional coverage argument is not claimed
novel, efficient, or fully justified for the fitted high-dimensional procedure.

## Candidate means and contrasts

For arm a, let Z_a,0 be the trusted target mean estimator and Z_a,j the
source-assisted mean from source j. Assume at least q_a of the K source means
are valid for mu_a; their identities may differ between arms. The target
counts as an additional valid candidate in each arm.

Form all (K+1)^2 contrasts

    D_jk = Z_1,j - Z_0,k,   j,k in {0,...,K}.

At least Q=(q_1+1)(q_0+1) contrasts are valid for TATE. This can use a source's
compatible control arm even when its treated arm is incompatible. It does
not require the same source to be valid in both arms.

Obtain the standard error of each contrast from the FULL joint covariance
of the arm-specific candidates:

    v_jk = Var(Z_1,j) + Var(Z_0,k) - 2*Cov(Z_1,j,Z_0,k).

The covariance includes shared target predictions and, when j=k, the local
two-arm residual covariance. Neither independent-pair formulas nor a common
target-factor equality should be imposed on finite fitted outputs.

## A dependence-robust reference

Choose r<=Q before inference and allocate alpha=alpha_A+alpha_S. Let I_A be a
target-only TATE interval with failure probability at most alpha_A. Construct
each pair interval I_jk with valid-pair failure probability at most

    p = alpha_S*(Q-r+1)/Q.

For exactly standardized Gaussian marginals this uses critical value
Phi^{-1}(1-p/2). Define S_r as the set of values contained in at least r pair
intervals, and report I_A intersect S_r, with the same recorded empty-set
fallback and optional hull as the scalar reference.

Fix population-valid arm subsets of sizes q_1+1 and q_0+1. Their Cartesian
product contains Q valid pairs. If fewer than r of those pairs cover TATE,
at least Q-r+1 fail. Markov's inequality bounds that event by
Q*p/(Q-r+1)=alpha_S. Union with anchor failure proves coverage at least
1-alpha. Pair dependence is unrestricted in this argument; the pairs do
NOT represent (K+1)^2 independent samples.

When q_1=q_0=0, use the ordinary target-only alpha interval directly: the
borrowing assumptions provide no additional valid source information.

## Growing-K reduction and limitations

If the anchor error bound has approximation error epsilon_A and each valid
pair error bound has an additional error at most epsilon_pair, the same proof
gives coverage at least

    1-alpha-epsilon_A-Q*epsilon_pair/(Q-r+1).

Choosing r/Q bounded below one avoids a quadratic K multiplier in this
coverage reduction. It does not establish the required uniform marginal
normal approximation or variance consistency. Those require the calibrated
score and nuisance-rate work already identified. Uniform pair remainders can
be bounded by sums of the two arm remainders; they do not need separate fits
for every pair.

The construction may be conservative. It is not asserted to attain the
Gaussian lower bounds or oracle length. Estimating q_a by counting apparent
agreement and treating it as known is unjustified. If an entire arm's sources
are known to be valid by assumption, pooling that arm is a separate useful
comparison; it must not be inferred silently from the simulated truth labels.

## Implementation route

1. Reconstruct the joint arm score covariance from saved fits and verify that
   D_00 reproduces the target TATE, and each D_jj reproduces its saved assisted
   TATE and empirical score variance.
2. Verify local cross-arm moment messages against retained score vectors,
   including uneven fold sizes and between-fold mean contributions.
3. Reuse the existing endpoint sweep for pair-interval voting. The explicit
   pair list has O(K^2) size; sorting costs O(K^2 log K). Avoid allocating its
   O(K^4) covariance matrix when only marginal pair variances are required.
4. Check exact Gaussian references with arm-specific invalidity, before fitted
   n=1000 experiments and any coverage claim.
5. A next-paper candidate-only fitter can omit the inner tasks used solely to
   learn aggregation weights. Its outer nuisance coefficients and candidate
   summaries must first match the current full fits. The production estimator
   continues to use its complete three-level weight-learning procedure.

The marginal version does not need the expensive variance-region Gaussian
integration or a common target-factor estimate. Such covariance-aware
refinements remain possible comparisons after their validity is established.

## Implemented evidence and remaining limitation

`arm_candidate_score_summary.R` reconstructs the full arm covariance and
checks local cross-arm summary pooling. `arm_pair_confidence_sets.R` forms
the pair intervals without allocating their K^2-by-K^2 covariance matrix.
The current covariance-input validation additionally uses an O(K^3) PSD
eigenvalue check on the 2(K+1)-dimensional arm matrix.

Twenty-six independent checks passed, including uneven folds, between-fold
cross-arm covariance, explicit contrast-matrix identities and zero-variance
Gaussian pairs. On 107 prespecified saved fits, maximum mean, covariance and
source-message identity errors were 1.11e-16, 1.58e-17 and 8.32e-20.

The first known-Gaussian study has 24 settings x1000 draws. Both pair quorum
choices had coverage96.6-98.2%. Partial-quorum mean interval lengths were
1.0458-1.1183 times the ordinary target interval length. This is a conservative
validity reference, not a satisfactory efficiency solution. Full quorum was
usually longer than target as well. Do not present these results as solving
the weak-deviation/growing-K efficiency problem.

`arm_pair_inference_reduction.md` gives the explicit error bound with estimated
SEs and nuisance remainders. Partial quorum avoids a K^2 multiplier in that
reduction; establishing its uniform fitted-score assumptions remains open.

# Structural weight candidate: smooth, rate-indexed quadratic bias penalty

This is a new candidate weight rule, not a claim about the current RoCE
implementation. It changes weight estimation, NOT an SE or confidence-level
coefficient. It is not adopted, and no manuscript or archived result changes.

## Candidate fixed before its sampling results are evaluated

With n equal to the target training sample size for that outer fold, minimize

    n * Var_hat(eta) + n^(3/4) * sum_j delta_hat_j^2 * eta_j^2.

The variance quadratic includes shared-target and cross-source covariance;
delta_hat is the training source-assisted minus target TATE discrepancy.
All treatment arms use the same TATE weights. The diagonal bias penalty does
not rely on cancellation between incompatible source biases. For fixed K,
the squared aggregate bias is bounded by K times the diagonal bias expression.

The exponent 3/4 is fixed as the midpoint of the sufficient rate window
(1/2,1) derived below, not chosen by RMSE or coverage. No exponent grid or
post-hoc cutoff search is authorized by this candidate protocol. The target
sample size anchors the objective scale; with fixed positive site proportions
it is asymptotically proportional to total sample size. Scaling the outcome
scales both variance and penalty by the same square, leaving weights invariant.

## Conditional rate argument, not yet a complete RoCE theorem

Suppose the variance quadratic has consistent O(1) coefficients after scaling
by n, its limiting curvature is positive definite, and there are finitely many
sources/folds. Suppose compatible-source discrepancies are O_p(n^-1/2), and
incompatible-source discrepancies converge to fixed nonzero values. The
inner discrepancies must describe the same limiting estimators actually used
on the outer fold; that is an explicit requirement, not automatic in v19.

For a general lambda_n, the compatible-source penalty is O_p(lambda_n/n), so
lambda_n/n -> 0 makes it vanish. Partitioning the linear normal equations
into compatible and incompatible coordinates gives incompatible weights of
order O_p(1/lambda_n), under the curvature/boundedness conditions. Their root-n
bias contribution is O_p(sqrt(n)/lambda_n). Thus the sufficient window is

    sqrt(n) << lambda_n << n.

lambda_n=n^(3/4) lies inside this window. Compatible weights converge to the
variance-minimizing compatible-source limit; incompatible weights need not
be exactly zero but their contribution is smaller than n^-1/2. Unlike the
fixed-cutoff Wald rule, null standardized noise does not retain an O(1)
random penalty in this limiting objective.

If the primitive estimators additionally have valid joint asymptotic linear
expansions and nuisance remainders, stable weights permit the ordinary
common-TATE influence expansion and its analytic variance. None of these
extra assumptions has been established for the full current p=100 pipeline.
This argument is pointwise under fixed separation, not uniform over local
alternatives. Finite-sample bias and coverage can still be poor.

## Implementation and review boundary

`quadratic_bias_weights.R` solves the small K-dimensional normal equations
using diagonal scaling and Cholesky solves. It never inserts an unreported
positive ridge or clips estimated weights. It also reconstructs the existing
common TATE pseudo-values from candidate fold weights using actual site sizes.

Only one candidate (power=3/4) may next be applied to all six rhos of the saved
C1/K2 n=100 diagnostic, retaining every seed and unchanged nuisance fits.
This is an exploratory reweighting comparison, not a new independent validation
set or formal production output. Report all settings and any failures; do not
change the power based on these results. Bootstrap is not run or substituted
for the analytic SE. The nuisance-score/target-anchor problems, inner/outer
alignment and full C1--C3/K2,4,8 experiment scope remain separate requirements.

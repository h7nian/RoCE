# Analytic weight-learning derivative: verified component, not full inference

The quadratic bias-weight candidate is unchanged (power=3/4). This work
derives a missing dependency in the candidate estimator, not a bootstrap or
a multiplicative SE adjustment. Nuisance fitted functions are still held
fixed; their derivatives/remainders remain a separate requirement.

## Derivative and observation bookkeeping

For each outer fold, write H eta + b = 0, with

    H = A + lambda_n diag(delta^2),  lambda_n = n_target_train^(3/4).

A is the target-training-n-scaled source-minus-target covariance matrix and
b is the corresponding target covariance vector. The analytic differential is

    d eta = -H^(-1) [d b + (d A) eta
                     + 2 lambda_n diag(delta * eta) d delta].

The implementation solves this small K-dimensional system; it does not form
a matrix inverse. For the estimator, the multiplier of d eta is the exact
global source increment L_k, containing target and source evaluation fold
fractions separately. Each observation's indirect derivative is accumulated
through EVERY training fold in which its original ID occurs. The direct
evaluation contribution is centered within its site and divided by that
site's full sample size. The indirect contribution is also site-centered by
the derivative identities, without treating overlapping appearances as new
independent observations.

Variance therefore decomposes as

    sum_i (g_direct,i + g_weight,i)^2
      = V_direct + V_weight + 2 sum_i g_direct,i g_weight,i.

The cross term is not dropped and may be negative. Within-inner-fold
centering of covariance components and globally pooled raw means follow the
actual aggregation code. Nominal site/training counts are fixed in this
empirical-mass functional. Deterministic positive-mass perturbations are used
only to verify its derivative, never to generate a resampling distribution.

## Real-data derivative checks

`check_quadratic_weight_influence.R` checks all six rho fits for seed10013.
All 48 deterministic directional checks pass; maximum error is 2.543e-12.
Ignoring the weight derivative yields an error as large as 0.00791935 in
these directions. Baseline moments, point estimates and fixed-weight variance
reconstruct exactly within the declared tolerances, and every site gradient
sums to zero. Outputs are in `quadratic_weight_derivative_seed10013_v1/` under
the independent-pilot root.

## Exploratory n=100 review

Job 18430464 completed 0:0 in 2:34. Every candidate point estimate and original
candidate fixed-weight variance matches the preceding reanalysis. No new
samples or nuisance fits were used. All output SHA checks pass under
`quadratic_weight_layer_n100_v1/`.

| Rho | Empirical SD | Mean fixed-weight SE | Mean weight-layer SE | Fixed-weight coverage | Weight-layer coverage |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.022648 | 0.022076 | 0.023146 | 0.93 | 0.93 |
| 0.5 | 0.025086 | 0.022476 | 0.025605 | 0.89 | 0.92 |
| 1 | 0.029027 | 0.023677 | 0.027893 | 0.83 | 0.93 |
| 1.5 | 0.029201 | 0.024530 | 0.028143 | 0.86 | 0.94 |
| 2 | 0.028904 | 0.024951 | 0.027991 | 0.89 | 0.94 |
| 2.5 | 0.028792 | 0.025197 | 0.027861 | 0.90 | 0.93 |

At rho1, mean indirect variance is 0.00007322 and the mean cross term is
0.00014740, so adding only a separate weight variance would miss most of
the change. At rho0 the cross term is negative in 44/100 seeds. This cannot
be represented honestly by a universal SE multiplier. The largest site
gradient-sum residual across 600 fits is 5.226e-17.

These are encouraging mechanism diagnostics on already known data, not a
passed independent statistical gate. Nuisance fits and fold allocation were
not differentiated by this n=100 calculation. In particular, ordinary mass
perturbations do not automatically encode the treatment-balanced fold design.

## First-order balanced-arm fold design

The actual partitioner allocates observations equally across folds separately
within A=0 and A=1, up to integer remainders. Conditional outcome/covariate
variation within each (fold,A) cell and the random whole-site treatment total
must be distinguished. If g_i is the complete layer's observation-mass
gradient and gbar_(a,k) is its cell mean, the first-order design map is

    B = (n_site / number_of_folds) * sum_k [gbar_(1,k)-gbar_(0,k)],
    g_design,i = g_i - gbar_(A_i,fold_i) + B*(A_i-p_site)/n_site.

The first term retains within-cell variation; the second represents the SAME
treatment-total fluctuation shared by all folds. Equal allocation per arm,
not independent fold-specific treatment fluctuations, determines B. Integer
remainder effects are not differentiated: their negligibility needs the
usual no-dominant-observation and regularity argument. This is a first-order
sampling-design calculation, not an exact finite-n distribution theorem.

`balanced_fold_gradient.R` has 12 passing assertions: a closed-form weighted
mean-A example, preservation of joint covariance/TATE contrasts, permutation
invariance and rejection of incompatible fold layouts. It was implemented
from the actual partition rule before reviewing the n=100 coverage table;
it is not selected by which SE gives better coverage.

The mapping also runs on every site/rho of seed10013, with actual arm/fold
count checks. For example, rho1's weight-layer SE changes from 0.0275795 to
0.0274249. The output is `quadratic_balanced_design_seed10013_v1/`. This is
not yet an all-seed design-mapped review, and nuisance derivatives are still
absent. Applying the map AFTER accumulating all paths preserves their cross
covariance; separately adding variance components would not suffice.

## Open requirements

The all-seed balanced-design review has now completed with coverage 91--94%
and unchanged candidate points. See `BALANCED_WEIGHT_LAYER_N100_REVIEW.md`,
including the less favorable one-case changes. Complete target CV states
have also been recovered with prediction parity, and the actual 603-dimensional
target system checked; see `TARGET_NUISANCE_SYSTEM_REVIEW.md`.

Do not promote the current table as full-estimator inference. Complete the
design-aware review, nuisance/target-anchor score construction, high-dimensional
rate and inner/outer alignment arguments, and independent validation before
the final simulation/RHC freeze. The original full experiment scope remains
unchanged. Original results and the root package are untouched.

Weight-gradient source SHA: `c05e5fbed603a737397271381213b938c5724635ac13828de98eba38d243691f`.
Design-map source SHA: `bcbace431b09d4636ece27f133a553c87a07fcdbb76d338e05b6d8e559d01c85`.
Design-map tests SHA: `99b287b4d52ac085a67794889ddfef33209e0e725506b6ff171d9b354396fec9`.

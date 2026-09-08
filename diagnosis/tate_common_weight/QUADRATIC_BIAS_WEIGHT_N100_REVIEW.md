# Quadratic bias-weight candidate: exploratory n=100 result

The candidate/protocol in `QUADRATIC_BIAS_WEIGHT_CANDIDATE.md` was fixed before
these reweighting results were examined. Job 18427523 completed 0:0 in 2:27,
using all 100 retained C1/K2 seeds and all six rhos. All 600 candidate attempts
succeeded. Original common-TATE point/variance reconstruction passed for every
fit before applying new weights. No nuisance model, clipping bound, observation,
bootstrap diagnostic or original output was changed.

Output:
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/quadratic_bias_weight_n100_v1/`.
Its five payload SHA checks pass. This is exploratory reuse of known data,
not independent confirmation or final production output.

## Point and interval results

| Rho | Candidate bias | Candidate RMSE | Original RMSE | Candidate coverage | Original coverage |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | -0.011478 | 0.025290 | 0.024953 | 0.93 | 0.93 |
| 0.5 | +0.006477 | 0.025786 | 0.025889 | 0.89 | 0.87 |
| 1 | +0.006886 | 0.029691 | 0.032112 | 0.83 | 0.82 |
| 1.5 | +0.004058 | 0.029336 | 0.033156 | 0.86 | 0.83 |
| 2 | +0.002217 | 0.028844 | 0.031728 | 0.89 | 0.90 |
| 2.5 | +0.001002 | 0.028665 | 0.030089 | 0.90 | 0.93 |

Target-only RMSE is 0.034767 and coverage 0.92. The candidate retains smaller
RMSE than target-only in every setting, but it is not uniformly better than
the original method: rho0 RMSE slightly increases, and coverage at rho2/2.5
decreases. At rho1/1.5, paired MSE differences versus the original are
-0.00014966 (MCSE 0.00003455) and -0.00023870 (MCSE 0.00004449), respectively.
These are exploratory paired diagnostics, not adjusted selection tests.

The interval calculation here deliberately remains the original realized-
weight pseudo-value variance with new common weights. It has NOT been
reinterpreted as the variance of the complete fitted candidate. At rho1,
empirical SD is 0.029027 versus mean SE 0.023677; at rho1.5 these are 0.029201
and 0.024530. Thus the candidate's conditional rate argument and improved point
performance do not, by themselves, close the finite-sample inference gap.

## Review decision

Do not promote this interval or declare the C1 coverage gate passed. Do not
change power=3/4 after seeing these outcomes. Continue derivation of the full
analytic effect of weight learning, and retain the separate nuisance-score,
target-anchor and inner/outer alignment requirements. No SE multiplier or
bootstrap replacement is authorized.

The candidate's smooth normal equations permit an analytic weight derivative:
if H eta + b = 0, then

    d eta = -H^(-1) [d b + (d H) eta].

This is a small K-dimensional linear solve, not a statistical rescaling.
The estimator derivative must combine its direct evaluation-score term with
the weight derivative from every overlapping training fold, indexed by the
same original observation IDs. The source increment in that chain rule uses
actual target/source fold fractions, not averaged fold contrasts. This
weight-layer derivative alone would still condition on nuisance fits; it
must not be advertised as the full estimator's influence function until the
remaining nuisance terms are justified or included analytically.

## Software evidence

All 20 candidate weight tests pass, including exact one-source solutions,
source permutation, shared covariance, scale invariance, rate examples and
invalid-curvature rejection. The unpenalized, zero-discrepancy K=2 solution
also matches installed v19's native optimizer within 7.566e-10 with unequal
source sizes. The main package remains unchanged.

Hashes at launch:

- Weight helper: `97f47e8ec7f747e115daed4193f22b88b717035c5bc709e070deb2d27e48df11`.
- Reanalysis: `26c3e186ea4deb4e986095f58ecd64f97caa9903bf6f75474ef8327981c61449`.
- Candidate protocol: `44820733c68218cdab0efee74b1e69e03ee97016f648b39b2ad9629b52d5ad94`.

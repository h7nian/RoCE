# Saved p=100 nuisance-system gate

This gate uses seed10013, rho0, outer fold1, source s1, treated arm from the
completed independent v19 pilot. It is a derivative/matrix/solver check, not
a new estimator, a selected regularization policy, or a coverage experiment.
No bootstrap or SE multiplier is used; no nuisance model is refitted.

## Actual system, not the simplified population fixture

`probe_saved_nuisance_projection.R` reconstructs the actual 2,010-dimensional
system: four initial target outcome fits, four initial source tilting fits,
one fold-summed final outcome fit and one fold-summed final tilting fit.
Each parameter block has 201 coordinates. It uses the original target/source
observation indices for the outer and inner folds, including unequal fold
sizes. Target/source loss normalizations are retained exactly.

The initial blocks use their out-of-two-fold training sets; the final losses
use all secondary calibration folds with their own initial plug-ins. The
final gamma target moment averages target fold summaries, whereas source
loss terms use the stacked source count. Target gradient summaries do not
apply M_fit truncation to the initial outcome logit; the source calibrated
loss does. The derivative implementation preserves that distinction.

## Completed checks

Slurm job 18424574 completed 0:0 in 20 seconds (1 CPU, 8 GiB requested).
Output: `results/direct_tate_mc500_b5000/independent_inference_pilot_v19/p100_projection_probe_v1/`.
All payload SHA checks pass. The driver verifies the frozen v19 installed
package and the exact seed bundle/artifact hashes before reading the states.

- Twelve deterministic directions check all 2,010 parameters jointly.
- Maximum analytic versus finite-difference Jacobian error: 5.095e-10.
- Maximum score-gradient error: 4.599e-11.
- Training score agrees exactly with native prediction/correction routines.
- Reconstructed final outcome/tilting KKT errors: 1.759e-7 / 3.286e-7.
- All three tested truncation counts in this unperturbed case are zero;
  nearest relevant boundary is 1.41878 away. Thus truncation cannot explain
  this case's finite fitted-score sensitivity.

Three solver fractions (1, 0.5, 0.1 of the maximum absolute training score
gradient) are fixed as numerical probes only, not chosen using truth or
coverage. All ten blocks solve at each fraction:

| Fraction | Explicit probe penalty | Max adjoint residual | Coefficient L1 norm |
| ---: | ---: | ---: | ---: |
| 1 | 0.065053281 | 0.065053281 | 0 |
| 0.5 | 0.032526640 | 0.032526649 | 1.39324 |
| 0.1 | 0.006505328 | 0.006505385 | 44.46110 |

Smaller penalty reduces the equation residual but produces much larger
coefficients. That is a reason to examine held-out function stability and
derive a rate-valid penalty policy, not proof of either failure or benefit.
Neither a production penalty nor a new CI is selected from this table.

## Active-truncation derivative branches

`check_saved_projection_truncation.R` uses copies of these fitted parameters
and adds six to one intercept at a time solely to activate the three relevant
branches. It does not pretend the modified coefficients solve a fitted model.
It then checks eight deterministic directions per branch and the native score:

| Branch | Activated source records | Max Jacobian error | Max score-gradient error |
| --- | ---: | ---: | ---: |
| Initial outcome plug-in | 201 | 5.477e-10 | 6.032e-11 |
| Initial tilting plug-in | 50 | 5.477e-10 | 6.032e-11 |
| Final inference tilting | 241 | 5.477e-10 | 6.120e-11 |

Native-score errors are zero. All original data, coefficients and published
results remain unchanged. These branch tests check the implemented piecewise
derivatives; they do not prove that truncation preserves identification or
asymptotic orthogonality under every data-generating distribution.

## Next implementation boundaries

The inner aggregation summaries use `per_k2_gamma` and `per_k2_alpha`, i.e.
initial out-of-two-fold fits, whereas the outer estimator uses the final
fold-summed calibrated fits. This is confirmed in
`R/cross_fitting_aggregation.R:420`. A score repair must handle both graphs;
applying the outer ten-block adjustment blindly to inner summaries would
misrepresent the fitted estimator and could use evaluation observations.

Still required: held-out score stability, both arms and target-anchor systems,
projection estimation/product-rate conditions, the federated summary contract,
and stable common TATE weights. This matrix gate alone does not resolve C1
coverage or authorize final production expansion.

The subsequent outer-fold evaluation is complete; see
`SAVED_PROJECTION_HOLDOUT_REVIEW.md`. The weak-penalty case passes numerical
equations but markedly increases held-out score variability. It is not
promoted, and the more moderate case is not selected using this holdout.

Hashes:

- Driver: `a8798c1244332fbf439476d23202babf6a26a997ae9a4f19700972dc48e7c0d5`.
- Solver: `7e149d700b3b9463fae7f21cb43c11a4b3c285568076b1795f92b2645a875999`.
- Truncation checker: `51068c276786712ce9f42d790ce208de0190129f9e01f0773db9df43d7061d27`.
- Original artifact: `9342f08bb428c1a6a62957b4110321b2b652710cb1b3a4c215490fd1da10c9ec`.

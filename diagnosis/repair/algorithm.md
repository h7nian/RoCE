# Three-level cross-fitting for the repaired RoCE method

The repaired method uses **three-level cross-fitting**. This refers to three
roles of the original fold partition, not three communication rounds. With ten
folds and 1000 observations per site:

| Role | Evaluation fold | Calibration folds | Initial-fit folds per block |
| --- | --- | --- | --- |
| Outer estimation | k1 | all except k1 (about 900/site) | exclude k1 and that calibration fold (about 800/site) |
| Inner validation | k2 | all except k1, k2 (about 800/site) | exclude k1, k2 and k3 (about 700/site) |

For a fixed training problem, each calibration fold contributes a loss using
its own initial nuisance predictions. One final OR vector and one final weight
vector minimize the combined losses. This happens separately for each site and
arm. Outer and inner source fits call the same source adapter; outer and inner
target fits call the same target adapter; both adapters use one loss assembler.

For each outer k1, all k2 validation means and within-fold second moments are
pooled with the existing sample-size normalization before computing Wald
discrepancies and optimizing weights. The final estimate is assembled from
site-aligned outer evaluation contributions, with each site's actual sample
size. No model coefficients are averaged across outer folds.

The three aggregation choices reuse these nuisance fits:

- `common_tate`: one source weight vector for the assisted TATE contrasts.
- `separate_arms`: separate objectives for the treated and control means.
- `joint_tate`: two weight vectors in one TATE-variance objective, including
  the covariance of both arms within each source and on the target.

The result records `crossfit_levels` independently of `communication_mode`.
The repaired candidate explicitly requests `target_nuisance_method="hou_calibrated"`
and `source_validation_method="calibrated"`. The package retains legacy
`lasso`/`initial` defaults for controlled comparisons; those executions are
labeled as two-level comparisons, not as the repaired method.

One-round source initialization uses target OR coefficients. Two-round source
initialization uses source OR coefficients, followed by the corresponding
target derivative moments. All fold-subset tasks can be enumerated before
weight learning. Verifying a complete summary-only inference protocol remains
necessary; naming a fit one-round is not a proof of that protocol.

The target PS loss uses the empirical moment
`mean((1 - I(A=a)) * initial_OR_derivative * augmented_features)`.
Its exponential term and the target OR loss are solved by the same kernels as
source calibration. The legacy recipe truncates initial OR predictors before
evaluating their derivatives and uses the initial tilt itself as the OR weight.
The experimental `score_derivative` recipe instead uses the derivative of the
untruncated OR prediction and the nonnegative negative derivative of the
truncated tilt. Its OR weight is zero outside the tilt truncation interval.
`calibration_control` also selects the initial target PS model and can specify
a target radius shared by fitting and evaluation. An explicit radius `log(9)`
matches the current DGP's [0.1, 0.9] target propensity bounds; it is not a
general-purpose recommendation. Source training and inference radii remain
separate arguments, and their statistical compatibility must be checked.
The stored empirical-score and active-set variance calculations do
not establish that clipping, nuisance remainders, or selection are negligible.
Those conditions and the Federated-DR/Pooled-DR baseline repairs still gate
formal claims about coverage.

The manuscript files in `docs/` describe the earlier estimator and its
unresolved assumptions. Their two-level formulas must be substantively revised
before representing this three-level procedure; a global terminology
replacement would incorrectly relabel the earlier derivation and Hou's own
two-level construction.

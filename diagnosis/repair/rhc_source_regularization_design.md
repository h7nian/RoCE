# RHC controlled source regularization study

The user authorized further work toward a stable RHC analysis. This exploratory
study follows the completed covariate comparisons, where representation changes
left the primary Medicaid weight ESS near25/193 and one record near60% of the
reported variance. The prespecified global1se alternative reduced concentration
but changed target and source fits together. Smaller SE or a significant effect
is not the criterion for selecting a new primary analysis.

## Fixed design

Use the historical61 working features and the same5039 records, Private target
and three insurance sources. Preserve confounders, outcome definitions, cohort
and imputation/scaling. Use Two-layer, one-round, joint TATE, cutoff2, radii5,
100 lambdas, proximal Newton, tolerance1e-10 and the verified CV certificate.
No data deletion, new overlap trimming, SE rescaling or DGP change is involved.
The ordinary and calibrated target references are reported separately.

Run all six variants under four outer partitions: the historical sample-size
partition and fold seeds101,202,303. This is24 model fits on one cohort, not24
independent patient samples or a known-truth coverage study.

| Variant | Global rule (target and initial OR messages) | Source initial weight | Source final weight | Source final OR |
|---|---|---|---|---|
| min | min | min | min | min |
| final_weight_1se | min | min | 1se | min |
| final_outcome_1se | min | min | min | 1se |
| initial_weight_1se | min | 1se | min | min |
| source_all_1se | min | 1se | 1se | 1se |
| global_1se | 1se | 1se | 1se | 1se |

`calibration_control$source_lambda_rules` controls the three source stages.
Omitted fields inherit the global nuisance rule. The target anchor and initial
OR messages continue using the global rule. For one-round estimation those
messages are target fits. The source controls also apply to the full nested
source calibration routine; they do not add communication rounds.

The final weight and OR losses use separate initial fits. Thus changing only
the final weight rule must leave source OR coefficients unchanged; changing
only the final OR rule must leave final weight coefficients unchanged. Changing
the initial source weight changes the final OR loss and the final-weight warm
start, and is interpreted as a coupled stage change. Initial OR rules are held fixed for all source-only
variants. No setting is selected separately for an influential patient/source.

## Validation and interpretation

Install a separate library copied from the checked production source snapshot;
leave the running simulation library untouched. The C++ kernels are unchanged.
Run control-validation, stage-isolation, nested-fold exclusion and installed-
package calibration/aggregation regressions. Before broad release, verify that
the new default reproduces the saved historical RHC primary fit. Compare both
target references and source model/prediction signatures within each partition.
Report and investigate any violation instead of silently accepting it.

Retain each result, stage rule, source coefficient signature, both-arm influence
functions, covariance, clipping and numerical diagnostics. Compare TATE and CIs,
source/arm ESS, largest weight and record variance share, ordinary covariate
balance, and variation across partitions. The same loss, training folds, CV grid
and tolerance apply; only the stated CV selection rules change. Exact numeric
lambda values can differ by site/arm/fold. Stronger regularization may trade
variance for bias; this observed cohort cannot establish lower MSE or coverage.

The completed four baseline fits on this same historical cohort remain the
reference comparisons; their fitting code and data are unchanged. Baselines are
not independently refitted for each RoCE tuning choice. The new main/default
fit must reproduce its matched ordinary-target benchmark before reuse is
accepted. Existing global1se and source-standard results remain named
sensitivities, not substitutes for this controlled comparison.

No manuscript result or global default is changed by this exploratory study.
Jobs and all new results are stored under the FACE-HD scratch RHC hierarchy.

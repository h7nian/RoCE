# Matched high-dimensional outcome CV result

Replacement array 18453192 completed frozen-v19 and corrected tasks with
exit code 0 in 1:03 and 1:13. Both used five native CV threads and returned
without warnings, invalid fold fits or final-fit nonconvergence. The two
`inputs.rds` files are byte-identical: same p100 outcome design, observations,
initial plug-in weights, 100-lambda grid and explicit CV fold IDs. All eight
result payload checksums passed.

| Package | Selected lambda.min | lambda.1se | Nonzero slopes in final fit |
| --- | ---: | ---: | ---: |
| Frozen v19 | 0.04353324 | 0.1329442 | 44 |
| Corrected CV training scale | 0.05754838 | 0.1329442 | 30 |

Independent minimum and 1SE recalculation from both complete CV paths agrees
with the stored selections; all candidates have five valid folds. The change
is consistent with removing the extra training fraction from the loss, not
with weakening regularization to obtain an attractive result.

The old shared-final C3 treated fit selected 0.1329442 under a different CV
randomization. The frozen package now selects 0.04353324 with the fixed new
labels. This is evidence of fold sensitivity in this input, not a claim that
the normalization correction alone resolves the earlier large penalty or
coverage problem. The two matched final fits have not yet been compared in
repeated-training TATE experiments.

Bundles `matched_v19_v1/` and `matched_corrected_v1/` are under
`results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/`. The original
parse-failure logs from 18453175 remain retained. No initial nuisance models,
TATE truth, source weights or SE multipliers were used to select these CV
penalties. Full estimator/variance integration and final experiments remain
outstanding.

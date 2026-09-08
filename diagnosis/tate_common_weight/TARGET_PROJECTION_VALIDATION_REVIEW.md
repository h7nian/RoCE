# Training-only target projection: bounded review

## Result and scope

Job 18436605 completed successfully using the procedure in
`TARGET_PROJECTION_VALIDATION_PROTOCOL.md`. This is seed 10013, outer fold 1
only, using cached target nuisance fits. It is method development, not an
independent coverage experiment or a production estimator change.

All four candidates completed all four inner validation folds. Recalculation
from the archived per-fold losses reproduced the candidate risks and selection:

| Projection penalty scale | Validation risk |
| --- | ---: |
| Null projection | 0 |
| 0.5 | 0.028143168 |
| 1 | -0.008880681 |
| 2 | 0 |

The selected penalty scale is 1. The archived identifier `c1` denotes this
scale, not simulation configuration C1. Selection used training-complement
validation risk, not outer-fold TATE, truth, coverage or interval width.
The outer-fold result was inspected after selection: coefficient L1 norm
0.2181998; target-fold TATE 0.1677125 to 0.1718610; score SD 0.9224372 to
0.9434397. These shifts do not establish improvement in bias or coverage.

## Reproducibility and code quality

The result bundle is
`results/direct_tate_mc500_b5000/target_nuisance_state_capture_v19/target_projection_validation_v1/`.
All six payload checksums passed. The sibling `target_projection_validation_v1_code/`
preserves source snapshots; all five entries in the result's `code_hashes.csv`
match these snapshots. The existing result bundle was not modified.

Subsequent helper cleanup names the argument `penalty_scale` explicitly,
validates nuisance dimensions and block partitions, and rejects invalid
derivative matrices and penalty inputs. The current driver additionally
reports `selected_penalty_scale`; the archived output retains its original
schema. The updated selector reproduces the archived selection and risks.
The target-validation suite passes 19 assertions and the shared sparse-solver
suite passes 31 assertions. These are bounded diagnostic tests, not a rerun of
the entire package test suite.

## Remaining gates

This target-only check does not repair the source ten-block nuisance graph,
establish the full fitted-estimator influence function, or validate final
TATE intervals. Extension must preserve exclusions for every outer and inner
fold, retain shared-arm and shared-target covariance, and undergo independent
statistical validation before production adoption. No bootstrap, SE multiplier,
production package change, manuscript edit or final simulation launch was made.

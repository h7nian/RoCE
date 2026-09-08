# Extend source validation to the remaining folds

Keep global outer fold 1, source s1, seed 10013 and both arms. The first
validation-fold-2 fit and equation audit passed. Extend validation to 3,4,5
without changing nuisance penalties, clipping or CV selection rules. Local
fold 1 maps to the chosen validation fold; local 2:4 map to the other three
global folds in increasing order. Global outer fold 1 remains absent.

The seed is 880102 + 10*(validation_fold-2) + arm. This preserves the first
split's seeds and avoids collisions across the eight fold/arm combinations.
The driver keeps validation fold 2 as its default. Original first-split code
is preserved in `source_validation_split_1_2_code_v1/` under the pilot root.

All four validation splits are required before selecting source projection
parameters. Fitting and equation checks alone do not select those parameters
or establish valid inference. No candidate may be selected from the first
successful fold alone, and failed remaining folds must not be omitted.

First-split audit 18439854 completed with exit code 0 in 16 seconds. Both
arms have 1,608 parameters; maximum Jacobian error 7.7271e-10, score-gradient
error 1.0604e-10, final KKT error 3.6175e-7 and held-out point error zero.
All three audit payload checksums passed. Array 18442505 was submitted for
validation folds 3--5 with at most two concurrent single-CPU tasks and a
40-minute operational time limit; no model parameter changed with this limit.

The generalized audit was rerun on validation fold 2 into
`source_validation_split_1_2_audit_v2/`. Both `equation_checks.csv` and
`validation_equation_states.rds` are byte-identical to their v1 counterparts,
verifying that fold parameterization preserves the first-split calculations.
Audit array 18442608 (tasks 3--5, concurrency 2) depends on successful
completion of fit array 18442505 and writes the corresponding per-fold audit
bundles. These later audits are not yet reported as passed.

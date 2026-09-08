# Complete outer-training C2/C3 source states

`fit_source_outer_basis.R` fits both source-s1 arms using global folds 2--5,
with fold 1 reserved for outer evaluation. Initial models exclude their
calibration fold as in the native one-round algorithm. Four calibration
pieces are combined by the frozen package's `process_source_site` entry point.
No selected projection coefficients are applied in this fit job.

Use seed 10013's fixed observations, rho 0, the verified working-basis
transformer, 100 nuisance lambdas, lambda.min, and M_fit=M_inference=5.
The new fixed fitting seeds are 881101+arm, shared across C2/C3 for the paired
comparison. These are new nuisance fits, not exact replay of the original
C1 fits. Do not replace a failed fit by its C1 counterpart.

Save global fold map 1:5 and `fit_scope=complete_outer_training`. The audit
supports this five-fold map separately from the four-fold validation views:
it must reconstruct four calibration pieces and each actual nuisance basis.
The resulting audited full-training states are required before evaluating
selected C2/C3 outer corrections. Stagewise hyperparameter selection must
still use the four separate validation complements, not these full-training
models on their own training observations.

The in-flight validation fitter is left unchanged. Its target-summary setup
and this outer fitter currently have parallel small adapters around the same
native calibration routine. Consolidate those adapters after the frozen
validation jobs finish, with identity tests; do not edit their running scripts.

Fit array 18444223 (tasks 2=C2, 3=C3) writes
`source_basis_C<config>_outer_1_v1/` under the independent-pilot root, where
the actual names are `source_basis_C2_outer_1_v1` and
`source_basis_C3_outer_1_v1`. Audit array 18444628 depends on its successful
completion and writes the matching `_outer_1_audit_v1` bundles. Submission
does not establish either fit or numerical audit success.

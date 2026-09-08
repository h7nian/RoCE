# Source projection validation risk: coupling constraint

The target nuisance Jacobian is negative block diagonal. The source nuisance
Jacobian is not: final calibration depends on initial nuisance fits. Therefore
the target's independent-block validation loss is not a joint source loss.

Write the source adjoint equation as J' a = D. For block i, hold the other
coefficients fixed and choose sign s_i=-1 for alpha blocks, +1 for gamma
blocks. The conditional convex quadratic uses

    H_i = s_i J_ii
    b_i = s_i (D_i - sum_{j != i} J_ji' a_j)
    L_i(a_i; a_-i) = 0.5 a_i' H_i a_i - b_i' a_i.

Its gradient is exactly s_i (J' a - D)_i. The existing backward block solve
has this interpretation. A conditional loss is not a claim that summing the
losses and jointly optimizing all coefficients solves the same equation:
the right-hand sides themselves depend on other coefficients.

`check_source_projection_risk.R` verified these identities on both saved
p=100, 2,010-dimensional source systems and all three probe settings. Across
60 block checks, conditional gradient error is at most 2.1e-16 and directional
error at most 9.2e-12. Conversely, the signed adjoint matrices S J' have maximum
asymmetry 0.9624048 (treated) and 0.5820288 (control). A naive joint quadratic
0.5 a' S J' a - (S D)' a differentiates the symmetric part, not S J' a - S D;
its gradient errors reach 0.0435445 and 0.0353611 on the nonzero probe fits.
The zero-coefficient candidate alone would not expose this mismatch.

The bundle is `p100_source_projection_risk_v1/` under
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
These are equation and finite-difference checks, not nuisance refits or a
selected validation policy. A stagewise held-out policy must fit and select
the downstream coefficients without using the current validation data before
scoring the corresponding conditional upstream losses. Alternatively, a
different justified joint criterion needs its own derivation. Neither simply
summing these losses nor dropping coupling terms is an approved shortcut.
The source split/refit graph and full analytic inference remain outstanding.

The first actual source validation refit is implemented in
`fit_source_validation_split.R` under `SOURCE_VALIDATION_FIT_PROTOCOL.md`.
Job 18439647 was submitted for both arms at outer 1/validation 2, writing
`source_validation_split_1_2_v1/` under the independent-pilot root. It reuses
the frozen native source-calibration entry point with restricted fold views;
no coefficient validation loss or production weight is selected by this job.
Inspect completion and numerical diagnostics before consuming its output.

The follow-on audit is `audit_source_validation_split.R`, submitted as
job 18439854 with `afterok:18439647`. It reconstructs the eight-block,
1,608-parameter system for each arm from the saved three-piece calibration
graph, verifies exact training/calibration/validation IDs, directional
derivatives, final-model KKT and held-out point identity. Intended output is
`source_validation_split_1_2_audit_v1/` under the same root. A submitted audit
is not a passed gate; check both task status and its committed output.

The fit job 18439647 completed with exit code 0 in 19:04 (extern cleanup
19:07). Both arms returned without failure or captured warnings. All four
fit-bundle payload checksums passed. At the subsequent status check audit
18439854 was still pending; final equation/KKT verification is not yet proven.

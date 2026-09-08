# Fresh analytic-weight pilot status

Seed 20001 completed in job 18457202 (exit 0, elapsed 02:24:37).
Independent numerical audit 18457617 completed (exit 0, elapsed 00:00:12).
The three simulation payload checksums and two audit payload checksums pass.
Results reside under
`results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/` in
`fresh_pilot_seed20001_v1/` and `fresh_pilot_seed20001_audit_v1/`.
The simulation manifest SHA256 is
`f95e48f01e74a05a9323f85083fdc2259ba4fb136c6c0f97d97a1d7261958b00`.

All six rho values and four method rows are present. Maximum algebra error
is 4.996004e-16; maximum weight derivative error is 3.222290e-12.
Selected-fit nonconvergence/degeneracy checks pass. The reported invalid CV
fold-fit count is 17091 at each rho; these reused diagnostics must not be
summed as six independent sets of failed fits. Rejected CV candidates remain
recorded, rather than being relabeled as universally successful fits.

This is one independent dataset, not a coverage or RMSE gate. For example,
at rho=1 the target-only estimate is 0.218453639, the original common-weight
estimate is 0.237431440, and the quadratic estimate is 0.226381078, against
truth 0.206265564. Both source-assisted estimates have larger absolute error
than target-only in this replicate. The quadratic fixed-weight SE is
0.022241869 and its analytic weight-layer SE is 0.029212364. These are
descriptive observations, not grounds for tuning or a full-influence-function
claim. Inference remains unvalidated and warning capture remains incomplete.

After inspecting this numerical gate, the remaining four prespecified seeds
were submitted unchanged, one replication per job, with dependent audits:

| Seed | Simulation job | After-success audit job |
| --- | --- | --- |
| 20002 | 18461760 | 18461764 |
| 20003 | 18461761 | 18461765 |
| 20004 | 18461762 | 18461766 |
| 20005 | 18461763 | 18461767 |

Poll these exact handles; do not duplicate submissions on observation timeouts.
Each output follows the same seed-specific directory convention as seed 20001.
Retain all six rhos and all failures. Five independent seeds cannot certify
nominal coverage. The full canonical simulation and RHC objective is unchanged.

# Control-arm source nuisance equations

The p=100 saved-state constructor now accepts an explicit `arm` argument,
defaulting to 1. The default constructed state was checked with `identical`
against the original treated-arm matrix-probe state. The two arm states use
identical target/source observation IDs and complementary arm indicators;
invalid arm arguments fail before construction. No production code changed.

Job 18438885 evaluated the control arm at seed 10013, rho 0, outer fold 1,
source s1. It completed with exit code 0 in 18 seconds. The retained bundle
is `p100_source_control_projection_probe_v1/` under
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
All five payload checksums passed.

The 2,010-parameter system includes all initial and final nuisance blocks.
Native score reconstruction error is zero; maximum directional Jacobian
error is 2.5902e-10 and score-gradient error 3.9602e-11. Native final-model
KKT checks passed. No inspected truncation branch is active at this fit.
All ten blocks solved for each of the three numerical probe penalties.
Coefficient L1 norms are 0, 3.792589 and 43.785255 for fractions 1, 0.5 and
0.1 of maximum score-gradient magnitude. These fractions are diagnostic,
not a selected source projection policy.

This establishes a control-arm equation check alongside the treated-arm
check, not a full TATE influence function. Their site-level scores must be
contrasted on the same observation IDs, retaining cross-arm covariance; the
two separate one-arm variances cannot simply be added. Source coefficient
selection, its inner estimation graph and joint target/source covariance
remain unresolved before final method adoption.

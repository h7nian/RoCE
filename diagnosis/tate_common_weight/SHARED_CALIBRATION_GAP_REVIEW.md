# Exact gradient-gap decomposition on calibration data

`decompose_shared_calibration_gap.R` decomposes the final score gradients
coordinate by coordinate. For gamma, with outcome residual r=Y-m_final,

    D_gamma = outcome_loss_gradient
              - mean_source[A U (w_final_clipped-w_initial_clipped) r]
              + mean_source[A U w_final_clipped r * clipping_active].

For alpha it separates the gamma calibration-loss residual, target and
source initial/final outcome-derivative differences, target fold-normalization
difference, and final inference-weight clipping. All eight arm/configuration/
parameter-block identities passed within 1e-12. Full coordinate vectors are
retained; maxima of individual components must not be added as signed effects.

Gamma-gradient component maximum absolute values:

| Config / arm | Outcome calibration residual | Initial/final weight gap | Clipping term |
| --- | ---: | ---: | ---: |
| C2 / control | 0.04391588 | 0.01768403 | 0 |
| C2 / treated | 0.05269628 | 0.02756716 | 0 |
| C3 / control | 0.04015512 | 0.04596941 | 0 |
| C3 / treated | 0.13294423 | 0.04775217 | 0 |

These are training/calibration quantities, not the independent population
gradients from the previous review. A small penalized KKT error does not imply
a small unpenalized score derivative: the loss residual can remain at its
L1 penalty threshold. Weight mismatch is a separate contribution. Before
another method change, inspect CV loss/penalty scaling and selected penalties;
do not infer a native solver error merely from these nonzero residuals.

The bundle `shared_calibration_gap_v1/` under the independent-pilot root has
three verified payload checksums. No fits, selected penalties or intervals
were changed. The shared-basis candidate remains unapproved for production.

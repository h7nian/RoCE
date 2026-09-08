# Selected source correction: residual sensitivity remains

The selected C2/C3 corrections were audited with coefficients held fixed.
For each nuisance graph the adjusted score derivative is D - J' a. All 48
directional checks agree with direct finite differences within 1e-7; 80
block-level sensitivities are retained in `selected_source_sensitivity_v1/`
under the independent-pilot root. All three payload checksums passed.

| Config / arm | Training max gradient, before / after | Outer evaluation, before / after |
| --- | --- | --- |
| C2 / control | 0.143924 / 0.057399 | 0.363917 / 0.394605 |
| C2 / treated | 0.123982 / 0.070909 | 0.243537 / 0.251495 |
| C3 / control | 0.055969 / 0.051769 | 0.202691 / 0.207875 |
| C3 / treated | 0.074826 / 0.062367 | 0.144408 / 0.155794 |

For C3, the selected null initial blocks do not make the coupled derivative
zero. Initial gamma-block maximum residuals reach 0.000855 (control) and
0.006932 (treated) on training data; initial alpha-block residuals reach
0.001301 for the treated arm. Original score derivatives for initial
parameters were zero before introducing final-parameter corrections.

This establishes remaining finite-sample sensitivity, not a proof that the
candidate is asymptotically inconsistent. Validation max-gradient norms are
noisy and are not the predeclared risk criterion. They must not be used to
retune the grid after seeing the outer fold. Conversely, passing L1 KKT bands
or selecting null blocks cannot be represented as proof of orthogonality.
The current procedure is not ready to justify full fitted-estimator Wald
intervals; nuisance/projection remainder control and independent sampling
validation remain substantive requirements.

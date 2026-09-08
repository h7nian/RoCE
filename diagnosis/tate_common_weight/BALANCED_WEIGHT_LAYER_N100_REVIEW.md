# Balanced-arm weight-layer review: all 100 saved C1/K2 seeds

The design map follows the actual treatment-balanced partitioner, not
favorable coverage. Candidate points and power=3/4 are unchanged. This still
fixes nuisance fits and does not differentiate integer fold remainders.

First launch 18431705 failed 127 before R due to an executable-path typo;
it created no analysis output and its log is retained. The shared fixed-module
launcher was then used. Job 18431903 completed 0:0 in 3:02, with unchanged
scientific script and seeds.

All 600 point/previous-variance identities and site/arm fold-count contracts
pass. Variance decomposes into within-cell variation plus a common site
treatment-total component. Output hashes pass under
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/quadratic_balanced_weight_n100_v1/`.

| Rho | Empirical SD | Mean design-aware SE | Fixed-weight coverage | Mass-gradient coverage | Design-aware coverage |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.022648 | 0.023015 | 0.93 | 0.93 | 0.93 |
| 0.5 | 0.025086 | 0.025454 | 0.89 | 0.92 | 0.91 |
| 1 | 0.029027 | 0.027732 | 0.83 | 0.93 | 0.92 |
| 1.5 | 0.029201 | 0.027983 | 0.86 | 0.94 | 0.94 |
| 2 | 0.028904 | 0.027832 | 0.89 | 0.94 | 0.94 |
| 2.5 | 0.028792 | 0.027703 | 0.90 | 0.93 | 0.93 |

Mean SE slightly decreases. Rho0.5 and rho1 each lose one covered observation
relative to the mass-gradient calculation; these less favorable values are
retained. Sampling design, not proximity to 95%, determines the derivation.

No bootstrap, SE inflation, new MC observations or nuisance refits were used
in this review. Nuisance estimation, inner/outer alignment, high-dimensional
conditions and independent validation remain open. The 91--94% range is not
a completed statistical gate.

Review SHA: `a903b0af5f2dc7c5578dd4be92f0a2557d60d27bd437f3abea1e4775fc17f06a`.

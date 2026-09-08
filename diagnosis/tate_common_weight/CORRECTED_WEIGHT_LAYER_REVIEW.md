# Analytic weight-layer candidate on the corrected-CV replay

The existing quadratic-bias common-weight candidate (fixed power 0.75) and
its analytic empirical-mass derivative were evaluated on the new full C1,
seed10013, rho0 fit. No nuisance model was refitted, penalty retuned or
bootstrap draw taken. Eight deterministic derivative checks passed with
maximum error 2.01e-12, and point reconstruction agreed within 1e-12.

| Quantity | Value |
| --- | ---: |
| Original common-Wald TATE | 0.1864364 |
| Quadratic common-weight candidate TATE | 0.1880036 |
| Candidate fixed-weight SE | 0.02170603 |
| Analytic weight-layer SE | 0.02179985 |
| Treatment-balanced-design weight-layer SE | 0.02165846 |

The indirect weight variance is 2.39703e-5 and its cross term with the direct
score is -1.98888e-5. Retaining the covariance is essential; adding the
indirect variance as if independent would overstate this layer's variance.
The treatment-balanced map can also lower the result and is not an inflation
factor chosen for coverage.

All four payload checksums in `corrected_weight_layer_v1/` under the outcome
CV candidate root passed. These are fixed-nuisance weight-layer derivatives,
not the full fitted-estimator influence function. They do not establish
nuisance remainder control or coverage. A fresh, declared independent
multi-replication evaluation is needed before considering this candidate for
the final method; all source/configuration and RHC requirements remain intact.

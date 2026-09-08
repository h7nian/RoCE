# Paired source-assisted TATE projection check

`audit_source_tate_projection.R` evaluated the treated/control source probes
on identical outer-fold-1 target and source IDs at seed 10013, rho 0, s1.
All training and calibration records were checked disjoint from evaluation.
The uncorrected difference reproduces the saved source-assisted TATE to
1e-12. This is a target-population effect assisted by source data, not a local
source ATE. All three numerical probe settings were retained.

| Fraction of each arm's maximum gradient | Original TATE | Corrected TATE | Corrected target score variance | Corrected source score variance |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.2135295 | 0.2135295 | 0.02679430 | 1.51242043 |
| 0.5 | 0.2135295 | 0.2131612 | 0.05450472 | 1.74108494 |
| 0.1 | 0.2135295 | 0.3504971 | 1.41728735 | 14.93454131 |

Fractions refer to the existing arm-specific diagnostic penalties, not a
selected common TATE projection policy. Similarity to the original estimate
does not justify choosing fraction 0.5. Fraction 0.1 produces substantial
instability in this check and cannot be promoted without resolving it.

For each site/setting, the paired variance identity

    Var(mu1_score - mu0_score)
      = Var(mu1_score) + Var(mu0_score) - 2 Cov(mu1_score, mu0_score)

was verified directly. For example, the uncorrected target score variance
is 0.02679430, not the independent-arm sum 0.03769815: the cross-arm term
is -0.01090384. These are held-out empirical score variances, not the full
estimator variance or a valid CI. No claim is made that the original v19
TATE implementation omitted this covariance.

The retained bundle is
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/p100_source_tate_projection_v1/`.
All four payload checksums passed. It includes paired site scores, all 12
site/setting/version covariance decompositions, and all three point checks.
No nuisance refits, observation removal, penalty selection, weight changes
or production launches were performed. Source training-only coefficient
selection and joint target/source/weight inference remain required.

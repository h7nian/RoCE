# Inner target correction: exact block decomposition

`decompose_nested_target_projection.R` reconstructs all 20 saved inner
projection scores, with 60 nuisance-block summaries and 4,000 observation
records (the same 1,000 target observations appear in four outer-training
contexts; these are not 4,000 independent observations). It verifies original
and adjusted score identity and

    Var(original + adjustment) - Var(original)
      = Var(adjustment) + 2 Cov(original, adjustment).

The output is `nested_target_decomposition_v1/` under
`results/direct_tate_mc500_b5000/target_nuisance_state_capture_v19/`.
All four payload checksums passed. No observations were removed and no
nuisance models were refitted.

## Largest mean shifts

| Outer/evaluation pair | Propensity mean contribution | Treated outcome | Control outcome | Total |
| --- | ---: | ---: | ---: | ---: |
| 2/4 | -0.00044436 | 0.05964562 | -0.00056769 | 0.05863356 |
| 4/2 | -0.00001597 | 0.05049207 | 0.00004807 | 0.05052417 |

Neither evaluation fold has active propensity or outcome clipping. The
treated-outcome projection has respectively 8 and 7 active coefficients,
with L1 norms 1.609653 and 1.548636. This localizes the observed large shifts
to the treated-outcome correction, not to active evaluation clipping or the
propensity correction. It does not establish that the treated-outcome
correction is biased or invalid in repeated sampling.

Adjustment variance is 0.07541038 / 0.09865067, whereas twice its covariance
with the original score is 0.3521968 / 0.3821320. Their sums, 0.4276072 /
0.4807827, exactly reproduce the score-variance increases. Thus treating the
correction as independent noise would miss the largest part of this change.
These are empirical score variances, not a complete fitted-estimator variance.

The largest individual adjustment accounts for 14.75% / 21.01% of centered
adjustment sum of squares. Keep all observations; neither trimming these
records nor reducing the chosen scale after viewing these shifts is part of
the fixed procedure. Before promotion, source and target corrections must
share a consistent estimation graph and the full analytic covariance must
include their overlap. The result does not justify production launch.

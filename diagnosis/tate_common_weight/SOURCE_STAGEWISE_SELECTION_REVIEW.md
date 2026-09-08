# First complete source stagewise selection

All four validation fits (global folds 2--5, outer fold 1, source s1, seed
10013, rho 0) and both arms passed equation audits. Across these eight
systems, maximum final KKT error is 4.8998e-7, Jacobian directional error
7.7271e-10, and held-out point reconstruction error 1.1103e-16. All inspected
fit/audit payload checksums passed; captured warning counts were zero.

Selection job 18443364 completed with exit code 0 in 46 seconds. The full
128-row candidate table contains no final or initial fit failure. Independent
recalculation from its CSV reproduces the archived two-stage selection.

| Candidate | Final-stage risk | Initial-stage risk within chosen final branch |
| --- | ---: | ---: |
| null | 0 | 0 |
| 0.5 | 0.057940175 | 0.00003286860 |
| 1 | -0.005736928 | -0.0000002835572 |
| 2 | 0 | 0 |

Both selected scales are 1. On refitting the projection equations using the
complete outer-training data, **all coefficients are exactly zero in both
arms**. This is distinct from the nonzero coefficients observed in smaller
validation training subsets. Consequently outer-fold source-assisted means
remain 0.4979431 (control) and 0.7114726 (treated), with unchanged TATE
approximately 0.2135295. Do not force a nonzero correction or reduce the
selected penalties in response. The result neither proves that nuisance
uncertainty is absent nor establishes valid full TATE intervals.

The bundle is `source_stagewise_selection_v1/` under
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
All five payload checksums passed. This is one development seed, one source,
one outer fold and rho 0. Other sources/rhos, the full inner source graph,
integration with target corrections and learned common weights, asymptotic
conditions and independent coverage experiments remain incomplete. No
production method, manuscript, or experiment manifest was changed.

# Proposed optional effect-heterogeneity sensitivity experiment

This is a prospective extension of the fixed default DGP, not a replacement
for completed or running experiments. No simulation has yet been submitted
with the proposed control.

The current bounded DGP uses outcome slopes beta=(.25,.20,.10,.05) in both
arms and intercepts -.5,.5. Its conditional TATE varies little: the population
SD is about .0082 in C1/C3. The common target contribution then induces only
about .000064 correlation between the two calibrated source candidates.
Consequently the existing simulations provide little evidence about stronger
shared-target dependence, even though the first-seed fitted target components
are not exactly equal.

## Proposed control and invariants

Add a nonnegative scalar `effect_heterogeneity`, default 0, with direction
d=(.8,-.3,.15,-.08). Use outcome slopes
`beta - effect_heterogeneity*d/2` in arm 0 and
`beta + effect_heterogeneity*d/2` in arm 1. The intercept difference remains 1. The planned values
are 0,.5,1,2. These values retain four nonzero slopes per outcome model.

The target and compatible sources share these conditional outcome means.
Existing source-outcome deviations remain additional intercept shifts.
Covariate and treatment-assignment laws stay as specified in bounded_dgp.md.
Thus bounded features, the sparse working bases and the correctness branches
remain available. Population truth changes with this optional control but
remains fixed across seeds, dimensions, K and source shifts within each cell.

The zero setting must reproduce the current DGP's data and truth exactly.
Record the nonzero setting in configuration, result metadata and experimental
identity; it must never reuse a different setting's results or checkpoint.
Do not alter the frozen installed library or current campaign configurations.

## Population evidence

`probe_effect_heterogeneity.R` evaluates calibrated population limits on the
four active coordinates. It uses fractional-response quasi-GLM fits only to
solve the population score equations; finite outcomes remain binary.
The 24/32 quadrature comparison differs by at most 5.95e-13. Candidate mean
identification and population gradients pass, and fitted predictor bounds
stay below the existing radius 12 (largest about 4.25).

| Control | C1/C3 target TATE | C1 conditional-TATE SD | C1 candidate correlation |
|---:|---:|---:|---:|
| 0 | .238608 | .00822 | .000064 |
| .5 | .236257 | .09896 | .00914 |
| 1 | .229634 | .19314 | .03463 |
| 2 | .208115 | .35719 | .11849 |

C2 has its own fixed truths and population projection; at control 2 its
candidate correlation is about .0846. Full results and checks are in scratch
`implementation/r7/heterogeneity_design_probe_v1/`. Population calculations
do not establish finite-sample coverage or high-dimensional convergence rates.

## Implementation and release sequence

1. Implement the optional control and verify exact zero-setting compatibility,
   unchanged assignment laws, parameter forwarding and fixed quadrature truth.
2. Check C1–C3 at p=10, n=1000/site, K=2 and controls 1/2 with matched baselines.
   Keep ten folds, three levels, cutoff 2 and the existing numerical settings.
3. Add identical-population controls with BOTH covariate_shift=0 and
   source_treatment_scale=0. This separates heterogeneous treatment effects
   from the estimand mismatch of pooling different site populations.
4. After execution and integrity checks, start a focused p=100 panel at
   50 repeats/cell; extend p=200 and source-outcome deviations after those
   results are available. Release gates must not require favorable coverage.

SS/IVW aggregate site effects and need additional conditions to target TATE.
Their implementation can be correct while their estimand differs under
covariate shift plus heterogeneous effects. Baseline reports must state this,
and retain the separate working-model assumptions of federated/pooled DR.

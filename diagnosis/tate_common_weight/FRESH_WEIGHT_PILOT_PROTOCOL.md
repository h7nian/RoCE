# Fresh corrected-CV analytic-weight pilot

Reserve new seeds 20001--20005, disjoint from the development pilot seeds
10001--10100 and canonical production seeds 1--500. Start seed20001 only;
inspect payload/numerical gates before adding the remaining four. Each seed
is one Slurm job covering C1/K2/p100, 1000 observations/site and all six
rhos {0,0.5,1,1.5,2,2.5}; use the existing declared rho-reuse mechanism.

Use the corrected CV library, original one-round nuisance algorithm, five
outer folds, 100 lambda values, min rule and M_fit=M_inference=5. Record the
original common-Wald estimator, target-only baseline, quadratic candidate
with fixed-weight variance and the same quadratic point estimate with its
analytic, treatment-balanced weight-layer variance. Power 0.75 is unchanged.
Do not introduce post-fit nuisance projections or shared-final-basis fits.

No comparison methods or weight bootstrap diagnostics are requested. The
simulation API requires a comparison-bootstrap count >=2 even when unused;
pass 2 only to satisfy validation and install an abort guard on the resampling
implementation. This is not a two-draw bootstrap. Retain complete nuisance,
common-weight, pseudo-value and gradient artifacts. Warning capture remains
explicitly incomplete across arm workers; do not promote this pilot workflow
to the final warning-provenance-complete production workflow.

First gates are numerical completeness and honest descriptive bias/RMSE/
coverage summaries with Monte Carlo uncertainty. Five repetitions cannot
certify nominal coverage. This pilot does not replace the original full
C1--C3/K2,4,8/500-repetition simulation and RHC objective. Candidate nuisance
remainder/variance conditions remain unresolved and inference_validated=FALSE
must remain until independent statistical evidence is adequate.

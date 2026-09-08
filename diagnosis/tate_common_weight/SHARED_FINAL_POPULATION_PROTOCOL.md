# Fresh fixed-fit evaluation of shared final calibration

Use 25 batches of 2,000 covariates per site, K=2 and p=100, with new seeds
928731+b. The earlier projection-development evaluation used 918731+b;
do not reuse it as independent evidence for this new candidate. Hold every
initial/final fit fixed. Compare original final calibration and shared final
calibration on the same new X values within each batch. Integrate A/Y
conditionally using the known DGP only for evaluation, never for fitting.

Record paired batch score means, target truth, shared-final alpha/gamma
gradient means, and original/candidate M=5 clipping counts. Reuse the batch
summary implementation to report integration MCSE, not estimator SE. The
generic CSV field `corrected_mean` refers here to the shared-final-calibration
candidate, not the earlier post-fit projection; retain that metadata label.

The primitive score helper passed ten coordinate finite-difference checks
with active truncation, covering target and source scores. Population
evaluation still does not replace repeated-training bias, coverage, remainder
or full TATE/common-weight integration gates. Do not select penalties or
alter M after inspecting this evaluation.

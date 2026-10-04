# RHC target rotation and bias-risk diagnosis

The user requested checking potential bias and other target sites. Real RHC
data do not identify estimation bias: neither an assisted estimate nor a
target-only estimate is known truth. Changing the target changes the target
population and estimand, so differences across targets must not be labeled bias.

## Frozen comparison

Retain the same 5,039 patients, historical 61 features, outcome (30-day death),
RHC treatment definition, imputation and pooled feature standardization. Rotate
each retained insurance group into the target role; use the other three groups
as sources. No patients or features are selected using outcomes or significance.

| Target | Sample size | Current evidence |
|---|---:|---|
| Private | 1,698 | Reuse completed discovery fits and four baselines |
| Medicare | 1,458 | New paired study |
| Private & Medicare | 1,236 | New paired study |
| Medicaid | 647 | New paired study; smaller target precision is expected to differ |

For each new target, run `min`, `source_all_1se` and `global_1se` on the historical
default partition and seeds101/202/303: 12 fits per target. Run the four unchanged
baseline methods once on the default partition: sample_size, inverse_variance,
federated_dr and pooled_dr. All method fits also return ordinary and calibrated
target-only estimates. Total new scientific jobs: 3 x (12 + 4) = 48.

Keep Two-layer, one-round, joint TATE, aggregation cutoff2, radii5, 100 lambdas,
proximal Newton, tolerance1e-10 and the checked CV certificate. Source-stage1se
retains min for target fits and initial OR messages. Baselines retain their
existing fitting rules and 5,000 bootstrap draws. Baseline uncertainty is not
converted to a cross-partition standard deviation. Do not average target effects
from different populations into a single estimate.

## Checks and interpretation

Confirm that rotating labels preserves every patient's X,A,Y, each insurance
group's sample size, the OR/weight basis, and feature centering/scaling. Preserve
the original Private-target data hash exactly. Compare source-only variants
with identical target stages within each partition. Verify saved influence
contributions reconstruct the reported variance, and verify the four new
baselines share the default-partition ordinary target reference.

Report point estimates and intervals; candidate/target-anchor discrepancies;
both-arm source weights and Wald discrepancies; source/arm ESS; individual
variance concentration; and ordinary covariate balance. Inspect whether target
covariate means lie outside each observed source-arm feature range. Those are
empirical support warnings, not proof of population nonoverlap or quantified
bias. Changing the target can reveal asymmetric support/transport problems;
it cannot alone distinguish treatment-effect heterogeneity, confounding,
misspecification, transport failure and sampling noise.

Do not select the target with the narrowest or most significant interval. All
four populations and all prescribed fits remain visible. No clinical causal
truth, coverage or RMSE is inferred from one cohort. Known-truth semi-synthetic
experiments would require an explicit generating model and would assess bias
under that model, not identify the bias of the observed RHC analysis.

Reuse the immutable checked v2 estimator library. Generalize only the existing
driver's target-site argument and the existing diagnostics/reporting labels;
no estimator kernel or loss is changed. Store new output under
`/scratch.global/zhan9381/FACE-HD/real_data/rhc/target_rotation_v1/`, organized by
target and method/baseline role. One fit per Slurm job; checkpoint/requeue;
preempt/public/saffo-2tb. Preserve old results and manuscript settings.

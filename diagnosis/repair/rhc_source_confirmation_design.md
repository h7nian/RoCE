# RHC source regularization: fresh partition confirmation

The completed six-rule, four-partition diagnosis is exploratory. Before these
additional fits run, fix three variants (`min`, `source_all_1se`, `global_1se`)
and twenty new outer-fold seeds 401 through 420, inclusive. This is 60 fits on
one existing cohort. The seeds are distinct from 0/101/202/303 in discovery.
Do not tune rules per site, arm, patient or observed confidence interval.

Retain historical 61 features and all 5,039 retained patients, the same source
and target definitions, Two-layer / one-round / joint TATE, cutoff 2, 100-lambda
grid, radii 5, proximal Newton and tolerance 1e-10. `source_all_1se` keeps the
target and initial OR messages at min. Omitted stage rules inherit the global
rule. The three configurations are paired within each new partition.

Report every prescribed result, including failures. Compare paired SE ratios,
estimate ranges, maximum individual variance contributions across all sites,
source/arm ESS, weight caps and ordinary covariate balance. Verify fixed-target
stage signatures and reconstruct variance from saved arm influence functions.
Review whether reductions persist across new partitions and whether worse
balance or remaining tail influence changes the interpretation. Smaller SE
and exclusion of zero are not sufficient selection criteria. This study has
no known causal truth or independent cohorts; it does not estimate coverage,
RMSE or establish less bias. Keep discovery and confirmation summaries separate.

Use the checked v2 source-rule library; its only additional estimator change
relative to discovery v1 forwards the same rules through the reused-source
refit path (not used by these full-fit RHC jobs). Installed-package regressions
and RHC interface checks must pass before submission. All output belongs in
`real_data/rhc/source_regularization_confirmation_v1/` under FACE-HD scratch.
One fit per Slurm job, 4 CPUs, 8 GB, 1 hour, checkpoint/requeue, with preempt,
public partitions and saffo-2tb. Preserve all old studies and baselines.

No global default or manuscript result is replaced by this confirmation. The
four already completed baselines remain available for the same historical
cohort; changing RoCE's CV rule does not require refitting unchanged baselines.

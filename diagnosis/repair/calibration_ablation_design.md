# Final score-calibration ablation

This comparison addresses the missing direct calibration evidence identified
on 2026-09-21. It changes the source nuisance program and target anchor in a
two-by-two design. It does not replace the default bounded DGP or main method.

| Source nuisance program | Target anchor | Experiment |
|---|---|---|
| Final score calibrated | Hou calibrated | Existing full method |
| Final score calibrated | Ordinary lasso | Existing `target_lasso` |
| Standard | Hou calibrated | New `standard_source` |
| Standard | Ordinary lasso | New `standard_source_target_lasso` |

The standard source program fits the initial merged-weight balancing loss
and an ordinary source outcome model on the full allowed training subset.
Its outcome loss has unit observation weights. Initial balancing remains;
the ablation removes the final score-calibration updates. Each inner
validation fold evaluates the complete selected program, trained without
the outer or validation fold. This differs from the old `initial_validation`
ablation, which evaluates fold-specific initial plug-ins.

The public control is `calibration_control$source_nuisance_method`, with
values `calibrated` (unchanged default) and `standard`. Complete inner
validation is named `source_validation_method="complete"`; `"calibrated"`
remains a compatible name. The standard program currently supports one-round
communication. Calibration on the target retains three fold roles; removing
calibration at both sites leaves outer evaluation and inner validation.

Keep 1,000 observations per site, two sources, ten original folds, bounded
C1–C3 data with no outcome departure, four active coefficients, cutoff 2,
100 nuisance lambdas, and proximal Newton tolerance 1e-10. Both programs use
the existing CV grid construction and selection rule, with matching final-fit
CV seeds on the same training subset. Numeric lambda values may differ
because the loss and its lambda maximum change. Aggregation, source residual
evaluation, covariance assembly, and held-out data are shared code paths.

The new plan has six campaigns: each missing factorial cell gets one
prespecified repeat per scenario at p=10, followed by 50 per scenario at p=100
and p=200. That is six execution checks and 600 high-dimensional repeats.
Every repeat is one Slurm job with checkpoints and requeue. First-seed gates
check execution and provenance; they never select favorable coverage or RMSE.
Existing paired baselines are reused after matching data hashes and ordinary
target references, without submitting duplicate baseline jobs.

Report bias, RMSE, empirical SD, mean SE, SE/SD ratio, coverage, interval
width, and paired squared-error differences with Monte Carlo uncertainty.
Keep calibrated target-only visible. Fifty repeats cannot precisely resolve
small coverage differences. A comparison without calibration is not covered
automatically by the calibrated score's single-model-misspecification proof;
its actual inference performance must be reported as an ablation result.

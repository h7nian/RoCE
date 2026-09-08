# Isolated outcome CV training-scale candidate

Source is copied from the frozen v19 source stage into
`results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/RoCE`.
Only the two outcome CV training-scale expressions and their comments change:
calibrated and refined branches use n_arm/n_source rather than
n_arm_training_fold/n_source. Validation observation-count weighting and
the final fitting objective remain unchanged. The root worktree package and
frozen v19 library/source are not patched.

The regression fixture uses centered +/- pairs and unequal arm/fold sizes.
It compares native CV scores with independent Gaussian closed forms or
one-dimensional binomial likelihood roots across both arms and both CV paths.
The test must pass all 16 assertions after candidate compilation. It is not
a coverage test or evidence that the corrected normalization improves bias.

Build/test job 18452740 uses `build_outcome_cv_scale_candidate.sh` and writes
an isolated library under `outcome_cv_scale_candidate_v1/lib`. Inspect its
actual status and logs before claiming success. Full installed R tests and
R CMD check are still required before any integration. Preserve all original
method results and do not relabel their CV normalization retroactively.

Build 18452740 completed with exit code 0 in 1:24; all 16 targeted assertions
passed. Full installed-library tests are submitted as 18452927, and isolated
R CMD build/check as 18452929. Repository tests and the new scale regression
were copied into the candidate source for package checking. These full gates
are not yet claimed to have passed. Installation emitted standard Eigen
compiler warnings in its retained log; review package-check status separately.

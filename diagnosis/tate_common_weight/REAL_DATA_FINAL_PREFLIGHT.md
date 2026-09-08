# Final RHC analysis: read-only preflight

The final real-data requirement is **not complete**. This preflight checks
available data and the rerun contract without fitting a new model, accessing
external accounts, or relabeling archived results as current evidence.

## Available data and retained cohort

`inst/extdata/rhc.csv` is locally available (2,320,412 bytes; 5,735 observations,
63 columns). Its SHA-256 is
`9ef4ab578be4b40ad5d97d3a7e08ffdc1f9f76aeeefee51b4996e4221556f8e8`.
The workspace and checked v15/v18/v19 installed copies agree. A v18 cohort/split
preflight reproduced 5,039 retained observations and 61 design columns:

| Role | Insurance stratum | n |
| --- | --- | ---: |
| Target | Private | 1,698 |
| Source 1 | Medicare | 1,458 |
| Source 2 | Private & Medicare | 1,236 |
| Source 3 | Medicaid | 647 |

The pre-specified support screen removes `Medicare & Medicaid` and
`No insurance`; it is not selected from effect estimates or confidence intervals.
In the RHC loader API, K=4 denotes total retained sites including the target.
The returned simulation-style source count is K=3 and n_sites=4; these must not
be confused when creating resource plans or reporting results.

## Locked analysis settings and numerical gate

The current driver is `scripts/slurm/run_rhc_direct_tate.R` with its matching
Slurm wrapper. Settings include target Private, outcome death30, site variable
ninsclas, the two-level support recode above, 10 folds, 100 nuisance penalties,
seed 42, M_tau=M_tau_inference=5, and common TATE aggregation cutoff c=1.

Fit the complete minimum-CV-loss nuisance rule first. Run the strict nuisance,
initial density-ratio support and pseudo-value audits. If the pre-specified
numerical gate fails, retain that attempt as a named sensitivity and rerun
**both target and source nuisance pipelines** with the one-standard-error rule.
Do not mix target-min with source-1SE and call it a full-1SE analysis; do not
choose the rule from a favorable interval or effect size. Historical min runs
hit a density-ratio coefficient boundary and line-search failures, so the
fallback is a genuine pending possibility, not a presumed passed analysis.

Comparison methods are sample-size, inverse-variance, federated DR and pooled
DR, with 5,000 paired cross-arm multiplier draws for their reported bootstrap
SEs. Those comparison draws do **not** account for RoCE adaptive-weight
uncertainty. The current RHC driver does not call the new weight-relearning
diagnostic; final RoCE inference policy must be settled explicitly.

## Existing outputs do not meet the final requirement

- `results/real_data/` and legacy forest outputs contain earlier arm-specific
  analyses, not the final common-weight TATE analysis.
- `rhc_supported_primary_provenance_v14_openmp` retains a numerically failed
  common-TATE min analysis.
- `rhc_supported_1se_v14_openmp_final` is a historical mixed nuisance-rule run.
- `rhc_full1se_v15_b5000` has corrected full-1SE fitting but an older package
  and cutoff provenance.
- The `production_v35_c1/rhc_primary` and `rhc_1se` results supply intermediate
  cutoff-1 reaggregation numbers, not final v18/variance-frozen evidence. The
  archived full-1SE audit does exist under
  `rhc_1se/audits/20260825_104652_job16921953/` and passed with zero failed
  checks; the matching min audit failed two checks. Both identify package
  `f99574520ff3c7fdbd276814edfd150b9f31f9832fd6ee792bf2bd4f9ef0ec80`.
  Their existence is not a final v18 audit, but must not be incorrectly
  described as missing historical validation. Old method labels and package
  provenance must not be relabeled away.

None of the inspected historical RHC results has the v18 package fingerprint
`27ebfb159e6e7f6db955be810e3ef89724741a6f02f68e83730a8f808e774b45`.
Matching a manuscript number is not a substitute for a final validated rerun.

## Software launch gate added before the final rerun

`run_rhc_direct_tate.sh` now requires `ROCE_PACKAGE_CHECK_GATE` and the installed
library's `audit_tests_passed.txt`. Before loading R or creating the analysis
output directory, it checks passed statuses and the current installed package,
package source, test suite, test-driver and check-driver fingerprints. Missing,
stale or duplicated required fields fail closed. This enforces the existing
software precondition; it does not freeze the inference policy or authorize
the final RHC analysis by itself.

The standalone `test_rhc_launcher_gates.R` passes 27 isolated launcher cases,
using real fingerprint routines and shell stubs instead of fitting models.
The cases include all eight missing/mismatched required keys, missing files,
duplicate/conflicting fields, changed installed/source/test bytes and exact
restoration. A separate read-only preflight against the real v19 library and
matching R CMD check gate also passed, stopping at the module-load boundary.
No RHC model was fitted and no analysis output directory was created by that
real-gate check. Slurm may create its own logs before any shell gate executes.

Current v19 package/source fingerprints remain
`2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7`
and `8e17115d87a5cfa0ba9f4051c4018989eb9c3b96067d5f9eeece792edd4b82ac`.
The independent simulation's nine-file workflow remains
`9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4`.
Only the RHC launcher workflow changes; its existing fingerprint includes the
launcher itself. Historical RHC artifacts retain their original provenance.

## Remaining final steps

1. Freeze the estimator and inference policy after simulation variance gates.
2. Build and verify the final installed package and R CMD check gates, then
   use new immutable output roots for the complete min/conditional-1SE workflow.
3. Run `audit_rhc_direct_tate`, `audit_rhc_initial_dr_support`, nuisance-detail
   extraction and the pseudo-value/comparison checks on the actual final fit.
4. Run the pre-specified cutoff sensitivity set {1, 1.5, 2, 2.5, 3}, with the
   primary cutoff represented exactly once and matching fit provenance.
5. Verify an exact reproducibility rerun using the corresponding comparison
   script. If adaptive-weight diagnostics are adopted, save them separately
   from the comparison-method bootstrap and state their inference limitations.
6. Render tables/figures and update the active manuscript copies only from
   matching, passed final outputs; preserve JASA and archived artifacts.

There is no observed truth for this real-data analysis, so it cannot itself
validate confidence-interval coverage.

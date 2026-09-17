# Current state of RoCE (Robust Federated Causal Estimation)

Last updated: 2026-09-17. Active iteration: HISTORY #0015 (nuisance solver and
n_folds redesign); #0014 (v2 recovery) is closed because production v2 was
cancelled on 2026-09-17 at n=500 with coverage 87-92% at K=4/8.
Historical decisions and superseded acceptance criteria remain in HISTORY.md.

## Solver and speed (2026-09-17, HISTORY #0015)

- Every L1-penalized nuisance fit now uses the proximal-Newton solver by
  default; `ROCE_NUISANCE_SOLVER=coordinate_descent` restores the original
  coordinate descent for paired audits against earlier production runs, and
  the production rows record the solver in `nuisance_solver`. Paired
  validation shows identical CV-selected lambdas, coefficients within 6.5e-6
  and 6-12x lower CV time; the user accepted it as result-preserving on
  2026-09-17.
- Decisions: K=8 runs use 5 CV threads (80 CPUs); the next design uses ten
  outer folds; both are being sized by `pilot_nfolds_K8_20260917`.

## Method and experiment contract

- Primary method A: `run_tate_crossfit(communication_mode = "one_round")`,
  common TATE source weights, `screening_rule = "soft_penalty"`, cutoff 1
  (`aggregation_lambda = 1`). A remains the provisional primary method.
- Rule B (`quadratic_bias`) is the pre-specified sensitivity estimator,
  computed from the same nuisance fits; it does not replace A.
- Reported `se` includes the delta-method weight-learning contribution.
  `se_fixed_weights` remains a diagnostic. Migration #0005 is complete.
- Both nuisance working bases are `[X - kappa, X^2]`. True-mechanism
  X-dagger mixing has strength 0.75: C2 misspecifies outcome, C3 propensity,
  C4 both. C1 is unchanged. Migration #0006 is complete.
- Five outer folds, `nlambda_init = 100`, nuisance rule `min`,
  `M_tau = M_tau_inference = 5`; comparison bootstrap draws = 5000.
- Negative transfer: C1-C3 x K = 2/4/8 (9 blocks), treated-arm shift.
  Shared shift: C1 x K = 2/4/8 (3 blocks), both-arm shift. Both families use
  binary outcomes, p = 100, 1000 observations per site, and rho =
  0/0.5/1/1.5/2/2.5, with 500 replicates per cell.
- Method specification: `docs/main.tex` and `docs/supplemental.tex`.

## Verified results and their limits

The v1 production artifacts have 25 replicates in every one of the 12
blocks (1350 negative-transfer task CSVs and 450 shared-shift task CSVs).
The library was reinstalled during v1; its current fingerprint 0d0de95d...
no longer matches the results and reuse gates (4cbbf8b0...). V1 is stopped,
retained at `results/direct_tate_mc500_b5000/production_20260908_v1/`, and
excluded from v2 aggregation. No v2 numerical results exist yet.

The 100-seed C1/K2 pilot must be compared using the same SE definition:

| rho | A weight-layer coverage | B weight-layer coverage |
|---|---|---|
| 0 | 0.93 | 0.93 |
| 0.5 | 0.91 | 0.92 |
| 1 | 0.88 | 0.93 |
| 1.5 | 0.94 | 0.94 |
| 2 | 0.93 | 0.94 |
| 2.5 | 0.96 | 0.93 |

Sources: `diagnosis/out/weight_layer/v1/replay_summary.csv` and
`diagnosis/out/quadratic_bias_rule/v1/rho_summary.csv`. Neither rule
uniformly dominates on coverage. These are historical pilot results,
not the new production estimates.

Target-only diagnostics at p=100, K=4, n_target=1000, 500 replicates:
C1 coverage 0.946, mean SE / empirical SD 1.010; C3 coverage 0.936,
ratio 1.070. The overlapping 25 seeds match v1 production estimates and
SEs to floating-point precision. C2 has no equivalent 500-replicate
confirmation yet. Files live under `results/direct_tate_mc500_b5000/`
`target_remainder_C{1,3}_20260910/`; their archived filenames start with
`c3_target_remainder_` but their config columns distinguish C1 and C3.
Increasing replicate count measures coverage more precisely; it does not
change the estimator or guarantee nominal coverage.

## Active work and acceptance policy

The user approved #0014 on 2026-09-11:

1. Repair rule propagation through reaggregation and rho reuse, protect
   installed audit libraries, propagate submission failures, validate
   diagnostic resume provenance, and synchronize C++ declarations.
2. Run regression tests and full package checks, then freeze a new v2
   source/workflow and library. Re-run A1-A5 on that version.
3. Use fresh `production_20260911_v2` and `Rlib_production_20260911_v2`
   directories under `results/direct_tate_mc500_b5000/`. Restart seeds
   1-500; do not combine v1 and v2 results.
4. Advance through 1/5/10/25/50/100/200/300/400/500. Implementation,
   completeness and provenance failures pause the affected block.
   Statistical warnings are recorded at 50/100/500 and do not stop the
   fixed 500-replicate experiment. The earlier coverage go/no-go criteria
   are superseded by this user-approved policy; they remain in HISTORY.
5. Produce per-family/config/K/rho summaries and A-primary figures with
   paired B and target-only comparisons and Monte Carlo uncertainty.

Current status: all A1-A5 gates PASSED on the frozen v2 package.
A1 382288: zero test failures; A2 382289: Status: OK; A3 382290:
both manifests passed. Package fingerprint:
ad51045358861e41a461c2079ed7fdbdd29f52b91205d232da303fe74324f456.
A4 382451/382452: 31/31 checks, verified output/sidecar checksums.
A5 negative transfer 382492/382494: 14 rows x 756 fields, zero exact
mismatches; A5 shared shift 382491/382493: likewise zero mismatches.
Both six-rho validation groups also match archived v1 core numeric values
for all 14 methods on seed 1; this is a regression check, not a coverage claim.

Production progress (observed 2026-09-12T10:42:52.258021-05:00). Verified counts below require
the actual per-block QC gate with matching frozen fingerprints. Submitted
rungs are not completed replications. Job IDs and current targets are in
`results/direct_tate_mc500_b5000/v2_progress.json` and `ladder_jobs/`.

| Family | Config | K | Verified replicates per rho | Submitted rung |
|---|---|---|---|---|
| negative_transfer | C1 | 2 | 200 | 300 |
| negative_transfer | C1 | 4 | 200 | 300 |
| negative_transfer | C1 | 8 | 10 | 25 |
| negative_transfer | C2 | 2 | 100 | 200 |
| negative_transfer | C2 | 4 | 100 | 200 |
| negative_transfer | C2 | 8 | 10 | 25 |
| negative_transfer | C3 | 2 | 100 | 200 |
| negative_transfer | C3 | 4 | 100 | 200 |
| negative_transfer | C3 | 8 | 10 | 25 |
| shared_shift | C1 | 2 | 50 | 100 |
| shared_shift | C1 | 4 | 50 | 100 |
| shared_shift | C1 | 8 | 5 | 10 |

Each block advances only after its own implementation QC, through the fixed
500-replicate ladder. Stage reports at n=50/100/500 include A/B/target metrics,
paired MSE and coverage differences, and PDFs. Final family aggregation and
FACE-style primary-A figures follow all n=500 QC/report jobs. Run-specific
scripts and SHA-256 inventories live beside the results; the source snapshot
and installed package are not changed during production.

The standalone weight-rule precision diagnostic
is a reserve tool, not a second production run. RHC follows the simulation
summary. No cutoff, DGP, nuisance or SE tuning is authorized by low coverage.

## Operating rules and deferred work

- Never reinstall into a frozen production library. New builds use fresh
  directories. Re-test an installed package with `ROCE_TEST_INSTALLED=1`.
- Source/workflow fingerprints remain strict. Finish edits before gates;
  changes during production require a new version and new gates.
- R environment: `module load R/4.2.2-gcc-8.2.0-vp7tyde`,
  `R_LIBS_USER=/users/0/zhan9381/Rlibs`; expensive checks run on Slurm.
- Keep the existing manifests, per-K resources and same-seed rho reuse.
- Inner calibrated-nuisance sensitivity, continuous C2/C4 validation,
  worker-warning capture, API renaming, legacy cleanup, and one-vs-two-round
  comparisons remain deferred. No new method is part of v2 recovery.
- Overleaf `JASA/` is the advisor's review copy and is not modified.

To resume: read HISTORY #0014, inspect the v2 gate files and job logs,
query Slurm, and continue the highest incomplete stage above. Update this
status from artifacts rather than prior conversational claims.

<!-- v2-incidents-start -->
Recovery status (2026-09-12T10:34:03.723618-05:00):
- C2/K4 and C2/K2 seed 26: both resolved through frozen canonical independent fits, full n=50 validation, exact-byte publication, and replacement production QC. Both resumed n=100.
- Original failed/cancelled jobs and replacement provenance remain recorded in `v2_incidents/` and `v2_recovery/`.
- Fourteen stage reports have passed independent numerical and PDF review; see `results/direct_tate_mc500_b5000/production_20260911_v2/stage_report_audits/README.md`.
- C2 data-only preflight 487252 is fully audited: all 1,500 planned config/K seed records and 7,500 positive-rho comparisons passed record/provenance checks. Ten groups cannot reuse (K2: 26/282; K4: 26/105/160/296/380; K8: 26/105/160). Two are already recovered and eight have registered canonical chains; no invalid case is unhandled. See `rho_preflight_v2/case_reconciliation.json`. This is separate from the 500-replication estimator experiment.
- C2/K4 n=200 continuation 481666 completed after verified publication of canonical seeds 105 and 160; group 498743 is running, with cumulative QC 498744 pending. Release job 493168 passed (6s); receipt hashes and current scheduler state were checked. Eight proactive cases are registered; all five K2/K4 cases are published and verified, with two K8 cases still awaiting publication. See `proactive_canonical_v2/registry.json`.
- C2/K8 n=50 continuation 492860 now requires both n=25 QC 492859 and canonical seed-26 publication 493133. The frozen submitter correctly starts at the first uncommitted seed and still audits the complete n=50 prefix.
- Proactive C2/K4 seeds 105 and 160 passed per-case validation and exact-byte publication; both source/production CSV bundles and commit hashes were reverified. Cumulative n=200 QC is still pending.
- C2/K2 seed 282 is validated, published and byte-verified; its future n=300 cumulative QC remains required. C2/K4 n=300 continuation 498745 additionally requires seed-296 publication 493139, preserving its n=200 QC dependency.
- First shared-shift n=50 report (C1/K4) passed independent numerical/PDF audit. High-rho A coverage 0.76/0.76/0.80 is retained as an interim statistical review signal; A and parameters remain fixed through all 500 replications. See the combined stage-report audit index.
- C2/K4 seed 296 is validated, published and byte-verified; the publication prerequisite for n=300 continuation 498745 is satisfied, while its original n=200 QC prerequisite remains required.
- C2/K4 seed 380 is validated, published and byte-verified. All preflight-confirmed K2/K4 exceptions now have canonical results; ahead-of-rung results still require their full cumulative QC.
- Shared-shift C1/K2 n=50 also passed numerical/PDF audit. High-rho A coverage 0.70/0.74/0.82 and RMSE findings are retained as interim review signals under the fixed 500-replication design.
- First K8 canonical case, C2/K8 seed 105, is validated, published and byte-verified (493134/493135). Only K8 seeds 26 and 160 remain pending; future cumulative QC remains required.
- OPEN incident: C2/K2 seed 122 (502045_1622) failed at rho=2.5 because no CV lambda was valid across all folds. No group outputs were committed. QC 502046 and n=300 continuation 502047 are paused. Frozen canonical task-5622 probe 512362 is staged separately; diagnosis is ongoing. The prior data-only reuse preflight was valid for this seed, so this is a distinct issue. See `v2_incidents/C2_K2_seed122.json`.
Scientific parameters and A remain frozen. The full 500-replicate goal is incomplete.
<!-- v2-incidents-end -->

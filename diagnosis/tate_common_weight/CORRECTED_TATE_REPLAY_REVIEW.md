# Corrected-CV full TATE replay: bounded result

Job 18453587 completed with exit code 0 in 1:42:49 using fixed seed10013
C1/rho0 data, original outer folds, both sources/arms and the corrected CV
package. All three output payload checksums passed. Point and site-centered
variance reconstruction errors are exactly zero, as is target-only change.

| Version | TATE estimate | Analytic SE | Target-only estimate |
| --- | ---: | ---: | ---: |
| Frozen v19 | 0.18698444 | 0.02165514 | 0.16706525 |
| Corrected CV scale | 0.18643637 | 0.02167815 | 0.16706525 |

Truth for this retained experiment is 0.20626556. This single replay is not
a bias, RMSE or coverage experiment. Common fold weights are unchanged:
the scale correction affects final source outcome fits, whereas original
inner weight learning uses the initial source plug-ins. Preserve that
distinction when reviewing estimator/weight alignment.

All selected nuisance fits converged and update/threshold ratios are at most
one. No degenerate outcome fits or R warning conditions were recorded. However
the CV paths are **not failure-free**: initial density CV reports 8,663 invalid
fold fits (8,363 tail skips); calibrated density CV 8,401 (8,101 tail skips);
calibrated outcome CV seven invalid fits (four tail skips). These counters
are exactly unchanged from the frozen v19 fit. They concern rejected CV
candidates, not nonconvergence of the selected final models; neither hide them
nor impose an unexamined zero-invalid-path criterion after fitting.

Evidence is `tate_replay_seed10013_rho0_v1/` under the outcome-CV candidate
root. Full warning provenance remains explicitly incomplete. No bootstrap,
experimental projection/shared basis, or SE factor was used. Adaptive-weight
inference and finite-sample nuisance bias remain unresolved, and canonical
simulations/RHC are still incomplete.

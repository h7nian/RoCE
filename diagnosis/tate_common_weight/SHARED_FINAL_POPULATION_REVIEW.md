# Shared final calibration: independent high-dimensional check fails to improve bias

Evaluation 18451963 and summary 18451970 completed with exit code 0. All
evaluation/summary payload checksums passed. This uses fresh seeds 928731+b,
25 batches and 50,000 covariates per evaluated site, with fixed fits and
conditional integration over A/Y. No parameters were retuned.

| Initial configuration | Original fixed-fit error (MCSE) | Shared-final error (MCSE) | Paired change (MCSE) |
| --- | --- | --- | --- |
| C2 | -0.0039259 (0.0009694) | -0.0111836 (0.0010418) | -0.0072577 (0.0010586) |
| C3 | -0.0097895 (0.0014958) | -0.0128420 (0.0009785) | -0.0030524 (0.0009436) |

The integrated target TATE is 0.2059419 (MCSE 0.0003053). All shared-final
source scores had zero M=5 clipping among the 50,000 evaluated X rows per
arm/configuration. The original C2 treated score clipped one row; the other
original configurations/arms clipped none. Thus active inference truncation
does not explain the candidate's larger error in this evaluation.

Shared-final maximum absolute population-gradient coordinates are all in
the final tilting block: -0.1179593, -0.1077088, -0.1043532 and -0.1244668
for C2 control/treated and C3 control/treated. Their integration MCSEs are
0.0023339, 0.0022845, 0.0020178 and 0.0016342. The finite high-dimensional
fits do not attain the low-dimensional zero-penalty population behavior.

The candidate is not approved for production. The deterministic span repair
does not constitute a finite-sample guarantee under separately penalized,
initial-plug-in calibration. Do not lower penalties or select a favorable
configuration after seeing these evaluation values. Investigate the remaining
initial-versus-final calibration discrepancy and regularization remainder
before adopting any alternative. The original package, manuscript and all
unfavorable results remain preserved.

Evidence: `shared_final_population_v1/` and
`shared_final_population_summary_v1/` under the independent-pilot root.
These are fixed-fit expectations and integration errors, not repeated-training
RMSE or confidence-interval coverage. The full user goal remains incomplete.

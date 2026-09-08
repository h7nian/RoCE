# Direct calibration repair candidate: shared final nuisance directions

The selected post-fit projections retain substantial fixed-fit population
gradients. An alternative structural candidate is to make both final nuisance
fits span the union of the supplied outcome and tilting directions while
retaining their distinct initial working models. This is a proposed method
change, not a production modification or a manuscript update.

`check_joint_basis_calibration.R` uses the frozen native calibration kernels
at zero penalty with inactive truncation on two exact finite-support examples:

| Case | Original sensitive derivative | Shared-final-basis derivative maxima |
| --- | ---: | --- |
| Correct quadratic outcome, reduced initial tilt | outcome 0.04865675 | outcome 9.71e-10; tilt 2.66e-10 |
| Correct tilt, misspecified outcome | tilt 0.03996943 | outcome 7.45e-10; tilt 7.59e-10 |

Both repaired population point errors are zero to printed precision. In the
second example the true logit includes a cubic term outside the shared
quadratic final basis: maximum fitted outcome error remains 0.04722592.
Thus that check does not obtain orthogonality by making both models correct.

Why this addresses the earlier span problem: the final tilting calibration
can balance every final outcome derivative direction; the weighted final
outcome score can balance every final tilting residual direction. If the
initial outcome is correct, final outcome correctness supplies matching
derivative weights. If the initial tilt is correct and representable in the
shared basis, both tilts agree at the population solution and the weighted
outcome equation controls the tilting derivative. These statements still
require identification, inactive clipping and appropriate population limits.

This does not prove high-dimensional penalized rates, robustness under active
asymmetric truncation, or full TATE coverage. A shared basis changes the
final working models and communication cost; it must not silently relabel
old C2/C3 methods or reuse their old coefficient dimensions. Initial bases,
final bases and mappings must be explicit and tested before any integration.
The existing frozen package and all earlier candidate results remain intact.

The bundle `shared_final_basis_population_v1/` under the independent-pilot
root has two verified payload checksums. No bootstrap, SE factor, data
trimming or fitted-sample tuning was used in this deterministic check.

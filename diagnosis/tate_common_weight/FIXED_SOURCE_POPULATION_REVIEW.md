# Independent fixed-fit population integration

Jobs 18448239 (evaluation) and 18448464 (summary) completed with exit code
0 in 30 and 5 seconds. Conditional averaging over A/Y passed its exact
identity check. All evaluation and summary payload checksums passed.
The evaluation uses 25 independent covariate batches, 50,000 X values at
each evaluated site, with all fitted/selected coefficients held fixed.

| Configuration | Original fixed-fit error (integration MCSE) | Corrected fixed-fit error (integration MCSE) | Paired change (MCSE) |
| --- | --- | --- | --- |
| C2 | -0.0038325 (0.0008584) | +0.0004885 (0.0011496) | +0.0043210 (0.0007561) |
| C3 | -0.0110501 (0.0009884) | -0.0062270 (0.0009948) | +0.0048231 (0.0001423) |

The integrated target TATE is 0.2057992 (MCSE 0.0002223). Corrected
source-assisted expectations are 0.2062877 (C2) and 0.1995723 (C3).
These are fixed-fit population expectations, not averages over repeated
training samples. C3 retains a non-negligible error relative to integration
noise. C2's small fixed-fit error does not establish valid Wald inference.

Maximum absolute corrected population-gradient coordinates are:

- C2/control: gamma_final, -0.1394202 (MCSE 0.0027268);
- C2/treated: gamma_final, -0.1337597 (MCSE 0.0020791);
- C3/control: alpha_final, -0.1160711 (MCSE 0.0037614);
- C3/treated: gamma_final, -0.0593190 (MCSE 0.0007396).

Thus the residual sensitivities at these fitted parameters are not explained
away solely by the small outer evaluation fold's noise. This does not prove
an asymptotic failure at a nuisance probability limit, but it rules out a
claim that the current finite-fit correction has empirically established
near-zero population derivatives. Corrected point bias and orthogonality
must be assessed separately. No grid, coefficient, clipping rule or interval
was changed based on this evaluation.

Evidence is under `fixed_source_population_v1/` and
`fixed_source_population_summary_v1/` in the independent-pilot root. The
remaining calibration tangent-span issue, remainder control, complete
target/source/weight integration and repeated-training coverage gates remain
unresolved; formal simulation and RHC completion is not established.

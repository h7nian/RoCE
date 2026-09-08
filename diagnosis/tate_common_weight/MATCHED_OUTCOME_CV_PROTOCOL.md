# Matched p100 CV normalization comparison

Use the C3 treated-arm outcome calibration input from the saved shared-final
candidate: seed10013, outer1, source s1, rho0, 100 lambda candidates, M=5.
Keep X/Y/A, initial tilting predictions, warm start and lambda grid identical.
Generate explicit arm-only CV fold labels once by rule seed994102 and the
package fold-count routine. Supply these same labels to both frozen v19 and
the corrected package; retain complete inputs, CV paths, warnings and fits.

This is a mechanism check, not exact replay of the old shared-fit CV, whose
fold labels were not saved. Earlier original/shared-final fits also used
different fitting seeds; their paired fixed-fit evaluation measures the
resulting candidates, not the isolated causal effect of basis enlargement.
This new comparison isolates normalization from fold randomness instead.

Use five allocated CPUs and opt into native CV-fold parallelism; keep BLAS
single-threaded. No data resampling, initial nuisance refit or selection by
TATE truth is introduced. Input-bundle identity must pass before comparing
the paths. Full TATE statistical consequences remain a later gate.

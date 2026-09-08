# Shared final calibration: p100 fit and numerical audit

Fit array 18449770 completed C2 in 29:27 and C3 in 29:54, both exit code 0.
Both arms returned without failure or captured warnings. All eight fit-bundle
payload checksums passed. Audit array 18449926 completed C2 in 12 seconds and
C3 in 9 seconds; all six audit payload checksums passed.

The audit independently reconstructed the two final calibration objectives
from the **original** initial coefficients. Original initial alpha/gamma lists
were verified unchanged, rather than treating mapped representations as new
initial fits. Across four arm/configuration fits, maximum KKT residual is
5.431e-7, loss directional error 5.862e-11 and native residual-score
reconstruction error 2.776e-17. No initial-outcome or initial-tilting M=5
truncations occurred in the audited calibration rows. Inference clipping
diagnostics are retained with the site scores.

Outer-fold source-assisted TATE is approximately 0.2047246 for the C2-initial-
basis candidate and 0.2355558 for the C3-initial-basis candidate. These are
single-fold values, not error, RMSE or coverage estimates. The final dictionaries
are shared, so these candidates must not be relabeled as unchanged original
C2/C3 estimators.

Bundles `shared_final_C2_v1`, `shared_final_C3_v1` and their `_audit_v1`
counterparts are under the independent-pilot root. Independent fixed-fit
population evaluation, active-truncation analysis, high-dimensional remainder
control and full TATE/common-weight integration remain required before any
production adoption. No manuscript, production package or final experiment
manifest was changed.

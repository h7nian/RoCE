# First p100 shared-final-basis calibration fit

Use seed10013, rho0, source s1 and outer fold1. Reuse the actual C2/C3
initial nuisance coefficients from their complete-outer-training fit bundles;
do not refit or enrich those initial models. Preserve original coefficients
and mapped plug-in representations separately. Fit the shared dictionary
mapping using target training rows only and verify prediction parity on source
calibration rows.

Final alpha and gamma use the declared common dictionary (C2: site basis;
C3: outcome basis), 100-lambda CV with lambda.min, M_fit=5, and fitting seeds
883101+arm. Reuse the native fold-summed block-design and fitting routines.
Each gamma calibration target summary uses the original outcome prediction
derivative multiplied by the shared dictionary; alpha calibration uses the
original initial tilting predictions via their mapped representation.

Retain original/mapped initial parameters, final fits, warnings, failures,
fold IDs and basis specification. These fits require independent objective,
KKT, truncation and score checks before any statistical interpretation.
No inference intervals, production package changes or initial-model changes
are part of this bounded candidate experiment.

# Explicit shared-basis mapping for the calibration candidate

`shared_calibration_basis.R` maps both supplied initial dictionaries into an
explicitly declared existing final dictionary. It uses only training design
matrices, not outcomes, true-model fields or the known TATE. C2 declares the
site dictionary; C3 declares the outcome dictionary. It does not automatically
construct arbitrary feature unions or silently discard directions.

The mapping includes the intercept, so raw versus centered linear features
are represented without changing initial linear predictors. Initial fits
remain in their original bases; mapped coefficients are plug-in prediction
representations, not new fitted parameters. This distinction must survive
provenance and later influence-function calculations.

Nine unit assertions pass: preserved predictors on fresh rows, either declared
final dictionary, rejection of insufficient spans, rank deficiency, automatic
basis guessing, changed held-out relationships and mismatched dimensions.
For the actual p=100 C2/C3 paired data, the maps were fitted using target
outer-training rows only, then verified on all three sites including held-out
rows. Both final dictionaries have 200 columns plus intercept.

Rank and span checks are explicit restrictions of this prototype. A setting
whose union is not contained in either supplied dictionary, or whose training
design lacks full column rank, requires another justified representation;
this helper must reject it rather than silently use a pseudo-inverse. No
production package, initial nuisance estimates, observations, or manuscript
was changed. Native high-dimensional final calibration and inference checks
remain outstanding for this new candidate.

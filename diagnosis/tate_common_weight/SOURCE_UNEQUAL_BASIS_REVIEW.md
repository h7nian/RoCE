# Unequal source nuisance dimensions: equation fixture check

The saved-state constructor no longer assumes both nuisance blocks have
201 coefficients. It checks outcome and site coefficient lengths against
their respective design matrices, and verifies target/source basis widths
agree within each nuisance family. This is diagnostic code, not a change to
the production basis or calibration method.

`check_source_unequal_bases.R` tests outcome/site dimensions 101/201 and
201/101 on both arms, reducing saved coefficient/design blocks consistently.
Across 24 directional checks, maximum Jacobian error is 4.28e-10 and score
gradient error 7.30e-11; native score reconstruction errors are zero. Both
payload checksums in `source_unequal_basis_fixture_v1/` under the pilot root
passed.

These fixtures are deliberately labeled as **not refitted C2/C3 experiments**.
They verify rectangular coupling dimensions, derivative signs and native
score alignment. They do not establish nuisance KKT conditions, double
robustness, calibration tangent-span adequacy or coverage under C2/C3.
In particular, the previously documented tangent-span counterexample is not
resolved merely by passing these dimension checks. Actual misspecified-basis
fits and full method alignment remain required.

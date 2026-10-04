# A sufficient KKT certificate for density CV

The seed28 p200 checkpoint profile found that about87% of recorded CV time was
in density/merged-weight fits. Every one of the1303 completed density fits in
that snapshot had invalid CV candidates and a skipped path tail. Final refits
took about7 seconds in total, versus over12000 recorded CV seconds. Those
durations span attempts and parallel arms and are not worker wall time; initial
glmnet fits did not have comparable timing attributes.

Keep the100-lambda path, existing CV tolerance and final refit unchanged. A
possible exact acceleration is to avoid a fit only when it cannot satisfy
the existing KKT acceptance condition. This is stronger than merely showing
that its mathematical objective has no finite minimizer: a loose numerical
tolerance could still accept a point of an unbounded objective.

For an intercept-augmented design X, nonnegative source weights w(gamma),
linear target moment t and an unpenalized intercept, the smooth score is

    g(gamma) = t - X' w(gamma).

The normalization and derivative weights are included in w. This form covers
initial/refined/calibrated density scores with positive clipped or unclipped
tilts. Let d be any direction satisfying Xd>=0. For any admissible L1
subgradient u (u_0=0, |u_j|<=lambda),

    d'(g+u) <= d't + lambda ||d_slopes||_1.

If the KKT infinity residual were at most epsilon, its left side could not be
below -epsilon ||d||_1. Therefore no candidate can meet that tolerance when

    lambda < (-d't - epsilon ||d||_1) / ||d_slopes||_1.

A cheap candidate direction comes from a previous iterate's slopes. Normalize
their L1 norm, then choose d_0=-min_i X_i,slopes d_slopes. This makes Xd>=0.
The prototype adds a conservative floating-point margin to that intercept and
subtracts another margin from the bound. An unusable direction supplies no
certificate. No LP solver or alteration of the penalty grid is required.

The prototype checks passed31 assertions, including multivariate separation
with no constant feature, interior support, rescaling, nonfinite directions,
and a nearly stationary point of an unbounded objective that must NOT be
rejected at the current tolerance. On one prespecified p200 synthetic CV
design, an iterate yielded a certified penalty bound.03349. At a penalty below
that bound, the existing proximal Newton routine failed to converge and an
independent gradient reproduced its KKT residual.

The validation compares entire descending CV paths with the certificate
on/off: all finite fold scores, valid/invalid masks, selected lambdas and final
coefficients must agree for both proximal Newton and coordinate descent.
A certified failure follows the existing
failure branch, including restoration of the previous warm start and the
same three-failure tail rule. A useful direction from a failed fit can be
retained for the certificate without using that fit as an optimization warm
start. CV tolerance is currently1e-4 even when the final refit uses1e-10;
the certificate must use the CV tolerance.

The argument defaults to the checked reference behavior during validation.
Its value enters nuisance-cache provenance: otherwise an
on/off timing comparison could silently read the other mode's cached fits.
The certificate is computed in native fold loops without invoking R from
OpenMP workers. Direction normalization and nonfinite-arithmetic rejection
have independent boundary checks.

No running worker or frozen experiment library was changed. The development
native CV functions now expose `use_kkt_certificate=FALSE`; final fits retain
their original solver. The shared wrapper only skips a subsequent CV fit when
a previously obtained direction certifies failure at that penalty. It restores
the original warm-start/failure behavior through the existing CV caller.
The native validation completed all24 cases, with full fold-score,
valid-mask, lambda-min/1se and final-coefficient agreement. The original job
timed out after23 cases; its dependent resume revalidated those records and
completed the last case with exit0. The checked native library also reproduced
the old4026-assertion regression summary. In the nine p200 cases per solver,
the summed reference/certificate CV time ratio was1.590 for proximal Newton
and1.925 for coordinate descent. These are descriptive CV timings for the
checked inputs, not end-to-end performance guarantees.

The R integration now uses `nuisance_cv_certificate=TRUE/FALSE/NULL`, scoped to
the enclosing fit or simulation call, with FALSE as the effective default.
Its effective value enters nuisance cache inputs and simulation provenance;
fold RNG seeds are unchanged. Refit/reaggregation metadata and socket/fork
workers retain it. The R integration v2 passed690 focused assertions and the
full4122-assertion regression suite. V1's dropped diagnostic attribute was
fixed in the existing lambda-support helper and its failed evidence retained.
The C2,p10 two-round ten-fold pipeline matches the previous library in all
aggregation modes and7280 stored nuisance records; it has no certified skips.
The full-grid C3,p200 one-round pipeline also passed: old/off/on have identical
coefficients in7280 stored nuisance records and matching estimates, SEs and
weights in all three aggregation modes. The enabled run recorded8896 certified
failures. Explicit two-worker PSOCK-on validation passed as well. New repeat
campaigns can request `--nuisance-cv-certificate`; the launcher requires a
checked capability marker and hashes it into provenance. Existing frozen
experiments still use their original code.

Evidence: scratch implementation/r9/seed28_checkpoint_profile_v1/,
r9/density_cv_certificate_probe_v2/, r10/cv_certificate_checks_v1/ and
r10/cv_certificate_integration_v2/.

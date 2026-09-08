# Candidate test scope and negative control

The targeted scale test distinguishes the implementations: frozen v19 fails
all eight native-score comparisons while passing eight optimizer-status
assertions; the corrected candidate passes all 16. The frozen failures are an
intentional negative control, not newly introduced solver failures.

Installed test job 18452927 completed with no failed assertions or test errors,
but only 1,731 repository assertions ran and 28 tests were skipped. Inspection
identified CRAN-gated TATE, rho-reuse and parallel/group-propagation tests among
the skips. This run is not treated as the full installed-test gate. Its result
bundle and warnings (including small-class warnings in smoke fixtures) remain
under `full_tests_v1/` in the isolated candidate root.

The verification driver now explicitly sets `ROCE_TEST_INSTALLED=1` and
`NOT_CRAN=true`. Job 18452976 writes a separate `full_tests_v2/` and requests
four CPUs for the parallel tests. R CMD check 18452929 remains a separate
standard package check; its completion cannot replace the expanded installed
test run. All frozen v19 and primary experiment outputs remain untouched.

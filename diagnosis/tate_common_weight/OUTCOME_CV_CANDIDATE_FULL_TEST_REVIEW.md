# Outcome CV scale candidate: verified software gates

R CMD check job 18452929 completed with exit code 0 and exact `Status: OK`.
The expanded installed-library job 18452976 completed in 2:08 with:

- 2,293 repository assertions passed, zero failed assertions or test errors;
- 16 new outcome-CV scale assertions passed;
- 18 explicitly gated diagnostic tests skipped in that suite invocation.

Core TATE aggregation, rho reuse and nuisance-group propagation test files
were inspected and had no skipped tests. The relevant opt-in density-ratio
CV scaling audit was then run separately against the candidate source and
forced installed library; all seven assertions passed. The remaining gated
diagnostics are not claimed as executed. Small-class warnings from smoke
fixtures remain in the logs and result objects, not suppressed from reporting.

The full result bundle `full_tests_v2/` in the isolated candidate root has
four verified payload checksums. Installed package fingerprint:
`6a715f3a41f2f42ef58bc08310d4ecca061524b0ab555d52a6221101a63a0a93`.
The first suite run without NOT_CRAN remains preserved as `full_tests_v1/`
and is not substituted for this expanded test evidence.

These gates establish a tested normalization correction candidate, not a
solution to TATE coverage. Frozen v19 and root production sources remain
unchanged. Before integration, preserve warning provenance, quantify the
correction's statistical effect under matched settings, and resolve the
remaining estimator/variance issues; canonical simulations and RHC are still
incomplete.

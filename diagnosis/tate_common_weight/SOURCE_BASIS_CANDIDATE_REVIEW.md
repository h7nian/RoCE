# Coupled candidates on actual C2/C3 basis refits

The candidate evaluator now accepts the working configuration explicitly and
uses the same basis transformer as the fit/audit path. On both actual paired
refits, all three scales and both arms completed successfully. The full
corrected validation mean agrees with score minus moment/coefficient product
to 1e-10. No parameters were selected from these first-split evaluations.

Approximate source-assisted TATE contrasts from the reported arm means:

| Working basis | Original | Scale 0.5 | Scale 1 | Scale 2 |
| --- | ---: | ---: | ---: | ---: |
| C2 | 0.2189525 | 0.1493873 | 0.2122686 | 0.2189525 |
| C3 | 0.1659879 | 0.1371709 | 0.1727996 | 0.1659879 |

These are fixed-sample, one-source, one-validation-fold diagnostics, not
sampling bias or RMSE estimates. The scale-0.5 changes are substantial and
must not be hidden by selecting a more attractive result after evaluation.
Initial-coupling terms are included; differing upstream conditional risks
still require staged selection with fixed downstream branches.

Bundles `source_basis_C2_candidates_v1/` and `source_basis_C3_candidates_v1/`
under the independent-pilot root retain all candidate coefficients, site
score vectors and conditional risks. All eight payload checksums passed.
Neither a full C2/C3 selection nor a robustness/coverage gate is complete.
No bootstrap, SE multiplier or observation removal was introduced.

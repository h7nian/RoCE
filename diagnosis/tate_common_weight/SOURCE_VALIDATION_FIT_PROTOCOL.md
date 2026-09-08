# First source validation nuisance fit

Fix seed 10013, rho 0, source s1, global outer fold 1 and validation fold 2.
Construct fold views mapping local 1:4 to global 2:5; global fold 1 is absent.
Call the frozen package's `process_source_site` with local main fold 1.
Its calibrated fits then use only global folds 3:5, with initial models on
the two-fold complements of each calibration fold. Neither global fold 1 nor
2 may enter any training/calibration input. Evaluation on global fold 2 is
performed only by the final source-score calculation, after fitting.

Reuse the one-round target-summary formulas and source fitting entry point,
including the within-call lambda cache and warm starts. Fit both arms with
100 lambdas, lambda.min, M_fit=M_inference=5, seeds 880102+arm. These are new
validation fits, not a claim to replay original saved fits. No projection
penalty is chosen in this run and no target-only cached model is substituted
for the source algorithm's initial outcome model.

Retain all coefficients, initial target summaries, fit diagnostics, warnings,
failures, and original training/calibration/evaluation IDs. A successful
function return is not sufficient: perform equation/KKT and exclusion audits
before using the fits for source projection validation. This first split
does not replace the other validation splits or full inner source graph.

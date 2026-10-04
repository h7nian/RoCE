# Fitted-candidate bias and covariance diagnostics

This research study addresses the gap between the checked population score
identities and a usable fitted-data confidence procedure. It does not alter
the current paper's DGP, estimator, experiment libraries or layer comparison.

The existing p4 quadrature audit is now a shared function. It checks each
specified source outcome shift against the generated conditional outcome law,
and retains nonzero private cross-arm covariance. Before new fits, its rho-zero
packets must reproduce the saved C1/C2/C3 packets exactly. Independent C1 oracle
calculations check treated-only and both-arm shifts, including the identity
Cov(I1 q1 residual1, I0 q0 residual0) = -E(residual1) E(residual0).

The initial fitted study is prespecified as follows:

| Control | Values |
|---|---|
| Observations | 1000 per site; 500 training and 500 independent evaluation |
| Covariates | p4, the four active bounded features |
| Scenario | C1, C2, C3 |
| Source count | K2, K8, K16, plus the target |
| Source deviation | rho0 or rho.5 in the first K/2 treated outcomes |
| Training seeds | 59301–59305 in every cell |
| Calibration | Three training blocks; one-round, score-derivative kernels |
| Optimization | Frozen checked library, proximal Newton, tolerance1e-10, grid100 |
| Candidate constructions | Translated source-specific OR and common target IPW OR |

This gives90 fitted repeats. The first C2,K8,rho.5,seed59301 case is an execution
and interpretation pilot; release is based on complete consistent outputs, not
small bias or favorable variance. Five training samples/cell are a diagnostic,
not a coverage experiment or proof of a growing-K rate. Varying finite K at
fixed n is likewise not a joint n,K asymptotic study.
Equal seed labels across K do not guarantee identical target observations;
the bounded DGP also changes source tilt positions with K. The paired
construction comparison is within each repeat, not between different K cells.

Each repeat saves fitted models, exact data/configuration hashes, checkpointed
nuisance fits, evaluation summaries and conditional population means/covariances.
Quadrature orders16 and24 must agree. The valid-coordinate labels come from the
simulation and are explicitly oracle information. Both candidate constructions
reuse the same fitted weights; common-outcome evaluation does not refit them.
Here common-outcome refers to sharing an outcome regression across sources.
The current paper's `common_tate` instead shares aggregation weights across
arms; this diagnostic uses oracle-valid GLS to assess fitted-score drift.

Report candidate mu1, mu0 and TATE conditional biases relative to their actual
conditional SDs; retain valid and invalid coordinates separately. Compare
estimated covariance with conditional population covariance. Oracle-valid GLS
summaries show the effect of knowing the compatible subset. When GLS uses an
evaluation-estimated covariance, the reported drift and variance hold its
realized coefficients fixed: they are not the conditional moments of that
data-dependent estimator. Do not silently label them as such.

This study does not supply a feasible nuisance-bias envelope, account for all
covariance-estimation uncertainty, or validate confidence intervals. Its purpose
is to determine which of those remaining errors matter in actual fitted data,
before claiming the Gaussian score reference extends to the fitted method.

One repeat is one Slurm job. Scratch-only outputs, completed-artifact hashes,
worker locks, atomic fitted checkpoints and requeue preserve reproducibility.
Current-paper experiments retain priority; start the pilot before releasing the
small research panel, with a modest concurrent workload.

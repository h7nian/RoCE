# C2/C3 full-validation stagewise selection

All six remaining paired validation fits and audits completed with exit code
0. Each new fit returned both arms without failure or captured warnings;
all inspected fit and audit payload checksums passed. Selection jobs
18444632_2 (C2) and _3 (C3) completed in 26 and 28 seconds.

Each configuration retained 128 candidate rows, with zero final/initial
candidate failures. Independent selection recalculation from each CSV agreed
with its archived selection state. All ten selection-bundle payload checksums
passed. The selected scales are:

| Configuration | Final calibration scale | Initial coupling scale |
| --- | ---: | --- |
| C2 | 0.5 | 1 |
| C3 | 1 | null correction |

C2 final-stage risks for null/0.5/1/2 are 0/-0.064112445/-0.008711476/0;
within its selected branch initial-stage risks are
0/0.0032797884/-0.0000164347/0. C3 final risks are
0/0.057956746/-0.006652543/0; initial risks within its selected branch are
0/0.00007815084/0.000001252926/0. The predeclared tie rule prefers null.

On the respective complete outer-training fits, C2 coefficient L1 norms are
3.232416 (control), 2.919267 (treated); C3 norms are 0.02212768 and 0.2683897.
The source-assisted outer-fold TATE changes approximately:

- C2: 0.2496970 to 0.2179182;
- C3: 0.2411156 to 0.2425866.

These are not RMSE or coverage results, and closeness to a desired TATE was
not used to choose the scales. C3's selected null initial correction is a
candidate choice, not proof that its initial-fit uncertainty is absent.
Its full adjoint residual and rate conditions still need statistical review.

Bundles `source_basis_C2_stagewise_v1/` and `source_basis_C3_stagewise_v1/`
are under `results/direct_tate_mc500_b5000/independent_inference_pilot_v19/`.
Scope remains seed 10013, rho 0, source s1, outer fold 1. Integration over all
sources, outer/inner folds and rhos, complete target/source/common-weight
covariance, independent coverage gates and the canonical simulation/RHC
experiments are not complete. No production package or manuscript changed.

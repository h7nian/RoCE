# Regularized nuisance-moment projection solver

Status: experimental structural-score component, not deployed RoCE inference.
No bootstrap, new MC sample, SE multiplier or positive statistical ridge is
used. No projection penalty has been selected for the actual p=100 method.

## Implemented objective and diagnostics

`sparse_moment_projection.R` solves

    min_a 0.5 a' H a - d' a + sum_j lambda_j |a_j|

by coordinate descent with cached gradients and a fresh gradient/KKT check
after every sweep. The penalty is supplied explicitly, either as a scalar
or per-coordinate vector. It is distinct from the aggregation penalty.
The solver never computes an inverse and allows positive semidefinite,
rank-deficient curvature when the penalized equations are supportable.

It validates matrix/vector dimensions, finiteness, symmetry and curvature;
normalizes a common objective scale; and reports KKT and original moment
residuals, objective, numerical rank and iteration count. Negative eigenvalues
within floating-point tolerance may be projected to zero with the adjustment
reported. No positive ridge is added to make an unsupported direction appear
identified. Original-system KKT is checked after any numerical adjustment.
Detected unsupported near-null directions and iteration-limit failures raise
typed errors with available diagnostic state, not successful coefficients.

`solve_joint_nuisance_projection` (internal dot-prefixed function) solves the
four-block adjoint system in reverse order, preserving the signs and coupling
to the initial outcome and tilting fits. The initial blocks are not silently
dropped. Nonzero regularization leaves bounded moment residuals rather than
exact orthogonality; controlling the corresponding statistical remainder is
a separate requirement.

## Verification completed

All 31 test assertions in `test_sparse_moment_projection.R` pass, with no
failures/errors/skips. They include diagonal closed forms, a dense-solve
zero-penalty reference, objective-scale invariance, a 60-coordinate/rank-20
problem, invalid and unsupported inputs, failure-state retention, and the
full coupled adjoint signs. The rank-deficient test checks the identified
fitted direction, not uniqueness of aliased coefficient representations.

The existing joint-observation TATE diagnostic now accepts `--sparse-equations`.
That path calls the blockwise solver six times (both arms, three model cases)
and verifies the entire adjusted score and site covariance identities at the
zero-penalty limit. Maximum TATE nuisance derivative is 8.327e-12; mean identity
error is 8.327e-17; covariance/site-scaling errors are at most 2.169e-19.
This is a low-dimensional integration check, not a positive-penalty p=100 gate.

Hashes at this checkpoint:

- Solver: `7e149d700b3b9463fae7f21cb43c11a4b3c285568076b1795f92b2645a875999`.
- Tests: `51c20bec92c518d9a1fdec58ab390be94a251c1b82eebaa9bb8d02ec8ea47416`.
- Nuisance fixture: `81fa23d9174690f788485dd84ec3a51285321cd01c1f54aee3fd1235b8dc9f58`.
- TATE integration: `31997b25e2fecdfcc3b353060f28f31181f6d954be265c949c22afc53f66ee88`.

## Actual saved-state mapping for the next gate

Read-only inspection of v19 seed10013/rho0 confirms that the needed fitted
states are retained at

    group_result$artifacts[[rho]]$direct_tate_results$one_round_crossfit$
      arm_results[[arm]]$fold_results[[outer_fold]]$source_results[[source]]

`gamma_s` and `alpha_ts` are the final calibrated vectors. Importantly,
`per_k2_gamma` and `per_k2_alpha` are the out-of-two-fold INITIAL vectors,
not separate final fits. `process_source_site` fits one fold-summed final
gamma and one fold-summed final alpha using fold-specific plug-in block
designs. The initial fits use complements of outer fold k1 and inner fold k2;
calibration observations come from k2. Do not replace that graph by an
unjustified arithmetic average of independently calibrated estimators.

The next real-data-matrix gate can reuse these states without rerunning all
nuisance CV. It must retain all per-k2 initial derivative blocks, the actual
fold/site normalizations, and active truncation derivatives. In particular,
the final gamma target moment averages the target fold summaries whereas
source loss contributions use the stacked-source sample size. Target initial
fits are shared across sources and must not be treated as independent.

Remaining requirements include a training-only projection-penalty policy and
high-dimensional rate/stability justification, the target-only anchor, the
federated summary/communication contract, and stable common TATE weights.
The C1 fixed-cutoff problem and full experiment scope remain unresolved.

The first actual p=100 saved-state gate is now complete (job 18424574, 0:0,
20 seconds). All 2,010 parameters, ten nuisance blocks, final native KKT
conditions and three numerical solver fractions were checked without refitting.
Active-truncation derivative branches were separately verified on coefficient
copies. See `SAVED_P100_PROJECTION_GATE.md` for exact results and limitations;
the large weak-penalty coefficient norm requires further held-out stability
review and is not a chosen inference policy.

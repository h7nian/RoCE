# High-dimensional validation and boundary experiments

Design decision: 2026-09-21, following the user's authorization to continue
the current paper's experiments and map their interpretation to the theory.
The plan was expanded after observing earlier deviation diagnostics; it is
not represented as a preregistered design. All results, including unfavorable
ones, remain in their original campaigns.

## Scientific questions

| Question | Experiment | Interpretation |
|---|---|---|
| Does the implementation behave as expected under the two nuisance-correctness branches? | C1 both correct; C2 outcome working model wrong; C3 assignment working model wrong; p100/p200 | Primary high-dimensional validation at the fixed site size |
| Does source aggregation reduce borrowing when the outcome mean is incompatible? | rho0/1/2, C1–C3 | Finite-sample assessment of separated, fixed source biases |
| What happens near the observed screening difficulty? | rho0.5, C1–C3 | Finite-sample boundary diagnostic; fixed rho is not a root-N local-alternative sequence |
| What remains when neither nuisance-correctness branch applies? | C4, rho0/1/2 | Explicit violation of the model double-robustness guarantee; no nominal-coverage promise |
| How sensitive is transport estimation to population composition? | covariate-shift multiplier0.5/2, rho0 | Change composition while keeping the outcome transportable |
| Can the basic pipeline behave sensibly without distributional or outcome differences? | C1, covariate-shift0, source-treatment-scale0, rho0 | Negative control for the transport component; both controls must be zero for identical covariate laws |
| What do the target calibration, inner calibrated validation and communication choice contribute? | Target Lasso; initial-only source validation; two-round protocol, each separately | Change one algorithmic component per ablation, retaining paired seeds |
| Does arm-specific borrowing or its optimization target matter? | common_tate, separate_arms, joint_tate within each method repeat | Share nuisance fits and generated data; retain cross-arm covariance |

The implementation and DGP remain frozen: n=1000 at every site, two sources,
sparsity4, bounded_joint_v3, ten folds, complete three-level calibration in the
main method, one-round main protocol, nuisance grid100, proximal Newton
tolerance1e-10 and the existing radii12. Initial-only source validation is an
intentional ablation; it is not a claim that the entire target/source pipeline
has only two levels. Current K remains fixed. Weak-deviation inference and
growing K are the next paper's direction, informed by Guo's JRSSB2018
TSHT/voting and JRSSB2023 searching/sampling papers.

The saved manuscript's oracle-selection theorem assumes separated biases
and a growing cutoff. The user's selected fixed cutoff2 is a finite-sample
configuration. Neither varying p at fixed n nor these cutoff2 results prove
that asymptotic theorem. C4 and deliberately biased source candidates have
their own interpretation and must not be presented as covered by the
pairwise identification theorem.

## Repeat budget and reuse

Primary C1–C3/rho0/1/2 cells at p100/p200 will each contain seeds1–200.
The rho0 main campaign is already complete. The original p100 deviation
campaign owns seeds1–50; continuation uses only seeds51–200. The new p200
deviation campaign owns seeds1–50 and includes rho0.5; its continuation owns
seeds51–200 for rho1/2 only.

Sensitivity, ablation and boundary cells use50 seeds initially. At nominal
95% coverage, the Monte Carlo SE is approximately1.54 percentage points for
200 repeats and3.08 points for50 repeats. Report uncertainty instead of
treating a point estimate near95% as proof of nominal coverage.

New work consists of4750 method repeats and2950 baseline repeats across20
campaigns. These counts include the paired negative control. The preparer
audits scientific experiment identities against every existing production
campaign before release. Different aggregation modes do not require separate
jobs. Baselines are shared across nuisance and communication ablations when
the underlying data, folds and baseline configuration are identical.

Baseline correctness and target-estimand validity are different checks.
SS/IVW average local effects and may differ from TATE under covariate shift;
the frozen federated/pooled comparators have their documented density-model
and inference conditions. Read `baseline_validation.md` when interpreting
those comparisons. The identical-population control helps distinguish these
estimand differences from an implementation problem.

All outputs belong under
`/scratch.global/zhan9381/FACE-HD/implementation/r7/highdim_validation_v2/`.
`profiles.json` defines the variations, `campaign_registry.json` locates old
and new blocks, and `design_checks.json` records the duplicate and baseline
coverage checks. No successful repeat is recomputed merely to obtain a
uniform directory layout.

## Release and analysis rules

Use the existing checked R library and repaired submission workflow. Each
repeat is a separate Slurm job with persistent exact-input fit checkpoints
and requeue. The sum of new campaign worker caps is1024 (768 method and256
baseline jobs); pending jobs count against each cap. The scheduler's shared
5000-job association limit is not treated as an entitlement to occupy every
slot. Existing and unrelated user jobs keep their identities and resources.

Each new campaign first runs one seed per cell. Gates require finite results,
correct data/target-reference identities and completed commits; they never
require favorable coverage, RMSE or weights. Baselines wait for the paired
method's initial execution checks. True numerical failures stop the affected
campaign for review. Known infrastructure failures are only retried after
review, with archived attempts and preserved seed/cache identities.

Analyze each scenario separately: bias and its Monte Carlo SE, RMSE,
coverage and its interval, empirical SD, mean SE, CI length and paired MSE
differences. Inspect treated/control source weights, discrepancies, penalty
activation, solver convergence and clipping diagnostics alongside these
metrics. Report method failures explicitly. Do not pool correlated scenarios
as independent repeats or choose new cutoffs after seeing favorable coverage.

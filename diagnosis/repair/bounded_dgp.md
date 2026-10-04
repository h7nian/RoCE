# Selected bounded simulation design

Version: `bounded_joint_v3`. This is the default DGP; historical `face` and
`roce` designs require explicit selection. The design was specified before
examining repeated-simulation coverage. No estimator is selected by these
initial diagnostics.

## Population and features

Each site has 1000 observations in the official validation profile. Target
coordinates are independent standard normals truncated to [-2.5, 2.5] and
standardized using their population variance. `p` is the actual number of
working covariates, with an additional intercept. There is no square expansion.
At most four slopes are nonzero in each true nuisance model. Intercepts and
combined nuisance supports still count separately in theoretical sparsity.

Target outcome logits have intercepts -0.5 and 0.5 and slopes
(0.25, 0.20, 0.10, 0.05). Target treatment slopes are
(0.35, -0.25, 0.125, -0.0625). Inactive coordinates are independent, centered,
and independent of treatment and active coordinates within every site.

For source j and arm a, let

    b[j,a] = delta * (j/K) * (0.25,-0.20,0.15,-0.10)
             + (2*a-1) * xi * theta / 2
    C[j] = E_target[exp(b[j,0]'u) + exp(b[j,1]'u)]
    f_source_j(x) * e_source_j,a(x) = f_target(x) * exp(b[j,a]'u(x)) / C[j]

Thus its ideal merged weight is `C[j] * exp(-b[j,a]'u(x))`.
Source covariates are sampled from this normalized mixture of tilted laws,
then actual treatment is drawn from its conditional propensity. The mixture
component is not used as actual treatment.

`dgp_control = list(covariate_shift = delta, source_treatment_scale = xi)`
has default (1,1). Delta ranges from 0 to 2; xi from 0 to 1. Identical source
and target covariate laws require delta=0 AND xi=0. Delta=0 with xi>0 retains
a small mixture-induced distribution difference.

## Scenarios and truth

C1 has correctly specified outcome and assignment models. C2 omits an
outcome interaction; C3 omits an assignment interaction. The wrong feature is
`u1 = X1*(1 + 0.75*X2)/sqrt(1 + 0.75^2)`, with the other three active features
unchanged. C4 is available as a supplementary both-wrong setting.
Source deviations can affect the treated arm or both arms, leaving the target
law unchanged. Changing source controls at fixed seed also preserves the
target's realized X, A and Y.

Population truth uses deterministic four-dimensional Gauss-Legendre
quadrature (24 nodes per coordinate), with 16-vs-24 convergence checks.
It does not change with seed, p, K, source shifts or source outcome deviations.
C1/C3 TATE is approximately 0.238608065665016. C2 has its own fixed truth.

## Relation to the theory and literature

The bounded sparse construction adapts the truncated-normal simulation in
[Tan (2020), Section 4 and Supplement I](https://arxiv.org/pdf/1801.09817).
Omitted-interaction scenarios also follow the experimental rationale in
[Hou, Mukherjee and Cai (2025), Section 5](https://jmlr.org/papers/volume26/23-1587/23-1587.pdf).
The explicit source-shift controls reflect the separate nuisance-correctness
and transport questions studied in
[FACE (Han et al., 2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12396575/).
The exact joint-tilt construction above is our adaptation, not a verbatim
simulation from any one paper.

Numerical population checks cover 960 active-block fits and 192 overlap
settings. Their largest fitted population predictor bound was 5.69666815;
M=6 and fitting/inference radius 12 are used in this validation profile.
These checks support the construction, not a complete inference proof.
The latest manuscript snapshot is revision
`8be3d007a28f6df9fb252a15f40d02b0bb083523`. Its revised inference theorem uses
a growing cutoff; the current simulation still uses the previously selected
fixed cutoff. Consequently DGP checks alone do not establish that theorem
for the current estimator, or establish 95% finite-sample coverage.

## Validation and reproducibility

Selected-profile single repeats for C1/C2/C3 at p10 all completed. The next
campaign is p10/20/50 x C1/C2/C3 x 200 seeds, K2, 1000/site, delta=xi=1, rho=0:
1800 independent jobs. Each uses ten folds, all three calibration/validation
roles, the full 100-lambda grid, proximal Newton at tolerance 1e-10, and exact
input caching. Common-TATE, separate-arm and joint-TATE weights share nuisance
fits and are compared against target-only estimation.

The controller defaults to 512 in-flight repeat jobs with 2 CPUs/job; pending
jobs count toward the cap. Its configurable safety cap is 2000, not the Slurm
account limit. First-seed checks cover all nine cells before further seeds.
Submissions are paced and batched independently of worker concurrency. No
failed or ambiguous submission is automatically retried. Full fit artifacts
are kept for each cell's first seed; result rows and failure records are kept
for every task. Reports include planned and successful repeats, bias, RMSE,
empirical SD, mean SE, coverage with Monte Carlo uncertainty, and paired MSE
comparisons. Low-dimensional results are reviewed before the p100/200 stage.

All outputs are below `/scratch.global/zhan9381/FACE-HD/implementation/r6/` or
the explicit campaign root. Frozen source, installed library and workflow
hashes identify the exact inputs. A scheduling change does not change seeds,
folds, estimator settings or convergence criteria.

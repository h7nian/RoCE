# Joint-nuisance orthogonal score: constructive prototype

This is an estimator-level candidate, not bootstrap, SE inflation, or a
truth-based bias adjustment. It is not adopted as RoCE and does not change
the original simulation, working-model definitions or manuscripts.

## Construction and signs

Stack all four fitted nuisance blocks in execution order:

    theta = (alpha_initial, gamma_initial, alpha_final, gamma_final).

Let M(theta)=E[m(O;theta)] be their unpenalized population estimating
equations. In the deterministic binary example these are target initial
outcome residuals, initial target/source tilting balance, source weighted
outcome residuals, and outcome-derivative-weighted tilting balance. The
last two equations depend on the initial fits; those derivative blocks are
included, not treated as fixed without justification.

For original source-assisted score psi and nuisance probability limit theta*,
write D = derivative E[psi(theta)] at theta* and J = derivative M(theta) at
theta*. If J is nonsingular, solve J' a = D and define

    psi_orth(O;theta,a) = psi(O;theta) - a' m(O;theta).

Then E[psi_orth(theta*,a*)]=E[psi(theta*)] because M(theta*)=0, while its
nuisance derivative is D-J'a*=0. Its derivative with respect to a is also
zero there. The adjustment is calculated from nuisance moment derivatives,
not chosen to match observed coverage or to increase an SE. It changes the
point-estimating equation. Its analytic variance must be calculated from
the complete adjusted score, not the old pseudo-values.

This follows the general influence-adjustment principle of
[locally robust semiparametric estimation](https://arxiv.org/abs/1608.00033).
High-dimensional estimation of the required representer needs additional
structure and rate control; a finite-dimensional matrix solve does not
establish those conditions ([Riesz regression reference](https://arxiv.org/abs/2104.14737)).
These general papers do not prove this proposed RoCE extension.

## Native-kernel population verification

`check_full_nuisance_orthogonalization.R` uses all four actual installed v19
native fitting kernels at zero penalty on exact expected Bernoulli losses.
Three support points make population expectations exactly evaluable. C2/C3
retain different outcome/tilting bases; no union-basis substitution is made.
The initial target treatment probability in this example is 1/2. Fractional
responses represent expected likelihoods, not fabricated binary simulations.
All truncation arguments are inside M=5.

| Example | Outcome / tilting dimensions | Original max nuisance derivative | Adjusted max derivative |
| --- | ---: | ---: | ---: |
| C1-style, both working models correct | 3 / 3 | 1.299e-9 | 5.551e-12 |
| C2-style, outcome model misspecified | 2 / 3 | 0.0600061 | 5.551e-12 |
| C3-style, tilting model misspecified | 3 / 2 | 0.0486568 | 5.551e-12 |

Analytic versus finite-difference Jacobians agree within 2.075e-11. Native
population moment residuals are below 5.737e-9. Projection equations close
within 2.082e-17 and all adjusted population means equal the target mean to
reported numerical precision. Thus the C2 counterpart of the earlier C3
tangent counterexample is also verified, and a full-nuisance analytic
correction removes these particular first-order sensitivities.

The earlier single-arm script revision had SHA-256
`382c58b7e57c4fac817cd0fdef5c63c633e6a805af067c9ef266fd9963c0f7c2`.
It has now been extended, without duplicating the fitting/Jacobian code, to
return both arm fixtures and allow treatment probability to depend on X.
Its default single-arm checks were rerun and reproduce the preceding results.
Current SHA-256: `77aa956c0589a4ff6870895c04ef7d07886dcb62e2baaa1f099744b4a5b2e755`.
The root package-source fingerprint remains `8e17115d87a5cfa0ba9f4051c4018989eb9c3b96067d5f9eeece792edd4b82ac`.

## Joint-observation TATE and covariance check

`check_joint_nuisance_tate_covariance.R` enumerates all 12 (X,A,Y) combinations
per site, using the same observations for both arms. Target treatment
probabilities are (0.3,0.5,0.7), and source treatment probabilities differ.
Both outcome means and treatment assignment depend on X. Source and target
covariate masses also differ. The calculation retains the C2/C3 feature-space
distinction and uses the full four-block correction separately for each arm.

For site-specific adjusted score contributions U_(a,r), the pairwise TATE is

    tau = sum_r E_r[U_(1,r)-U_(0,r)],
    Var_hat-scale = sum_r Var_r[U_(1,r)-U_(0,r)] / n_r
                  = sum_r {Var_r[U_(1,r)] + Var_r[U_(0,r)]
                           - 2 Cov_r[U_(1,r),U_(0,r)]} / n_r.

Using n_target=1000 and n_source=750 deliberately checks unequal site scaling.
The N/n_r-scaled, within-site-centered pseudo-value calculation gives the same
variance. The population TATE is 0.170903887208 in all three examples; maximum
mean identity error is 5.552e-17, full stacked TATE nuisance derivative is below
4.164e-12, arm covariance identity error is zero, and site-scaling error is
below 2.169e-19. Analytic variance is approximately 0.00131572237. Incorrectly
adding independent-arm variances gives approximately 0.00137764171 instead.
This is a check of the proposed score's covariance, NOT evidence that the old
v19 common-weight implementation omitted arm covariance.

Source-arm covariance is zero in these population examples because the
source residual contributions are arm-exclusive with zero source means.
That is an example property, not a general license to discard source-arm
covariance for finite fitted scores. The nonzero target-arm covariance is
retained. Likewise, nearly equal variances across these three small-support
examples are not a claim of general efficiency equivalence.

Both scripts execute successfully. TATE script SHA-256:
`553c3249fc601eaaeb409969826c5f6c14857414b0a144101f32690258899242`.
No data resampling, simulated replications, CI inflation, or production edits
were used. These are exact low-dimensional population calculations.

## What this does not yet prove

The prototype now also has a tested blockwise coordinate-descent projection
path (`--sparse-equations`) instead of a dense inverse/solve. Its zero-penalty
limit passes the same complete joint-observation TATE checks. See
`SPARSE_PROJECTION_IMPLEMENTATION.md` for 31 solver assertions, current file
hashes, the precise saved-state mapping and remaining high-dimensional gates.
The earlier hashes above describe preceding prototype revisions, not the
current solver-integrated files.

- These are finite-dimensional population examples, not the actual
  C1/C2/C3 p=100 DGP or a repeated-sampling gate.
- Penalized high-dimensional fits do not solve their sample moments exactly.
  Estimating a requires a justified regularized representer, stability bounds,
  and nuisance/projection product rates; directly inverting an empirical large
  Jacobian is not an accepted production implementation.
- The target-only anchor needs its own complete moment system and adjustment.
- The joint-observation check covers one target/source pair, not multiple
  sources with learned common weights. It does not verify the whole aggregated
  construction or a repeated-sampling TATE CI.
- The actual clipping/truncation operators and their derivatives must be
  included where active, or their contribution must be proven negligible.
- A federated implementation must declare the training/evaluation summaries
  required for coefficient estimation, mean aggregation and covariance, without
  silently adding a second communication round or ignoring communication cost.
- This correction does not remove true incompatible-source bias. The fixed-
  cutoff weight-stability problem remains separate and unresolved. A structural
  aggregation repair, not an SE multiplier, is still required for that problem.

No simulation or resampling jobs were submitted for this prototype. The full
primary/sensitivity/RHC requirements in COMPLETION_CHECKLIST.md remain open.

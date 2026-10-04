# Status of the calibration-evidence review

Updated 2026-09-23. Current Two-layer MC200 evidence is separated below from
the historical Three-layer factorial comparison and saved manuscript audit.

## Current Two-layer MC200 evidence

The p100 main source ablation, main target ablation and stronger-covariate-shift
source ablation each have600 verified pairs (C1-C3,200 seeds/cell,K2,rho0).
The complete overview is scratch
`implementation/r11/current_paper_mc200_v1/reviews/calibration_p100_overview_mc200_v1/`.
All data hashes, ordinary-target references and required component invariants
match. Both arm variances and their covariance reconstruct TATE variance.

Source calibration changes TATE RMSE by-0.24%/-0.58%/-0.03% in the main C1/C2/C3
cells, and-0.54%/-0.52%/+0.30% under covariate_shift2. Every total TATE paired
MSE contrast is within two MCSEs. Target calibration changes TATE RMSE by less
than0.3%, while reported/empirical variance ratios move from.848/.856/.934 to
.880/.893/.966. The paired changes in variance discrepancies are about3.0-3.9
jackknife SEs, but both versions still underestimate variance in these samples.
This supports neither a large isolated calibration gain nor uniform dominance.

The saved-score mechanism comparison is also complete:600 pairs at each shift,
2400 fitted programs reconstructed within1e-12 for both means, TATE and both
variance formulas. Actual outer folds agree; eta scores exclude held-out rows;
target anchors are invariant. Candidate and eta changes are small and sometimes
oppose each other. Both algebraic MSE decomposition paths check exactly. These
weight exchanges are diagnostics using both programs, not new estimators.

Five focused Python tests and three existing review tests passed. No estimator
or scientific input was changed. p200 continues; do not mix partial p200 output
or historical Three-layer results into these current paired MC200 conclusions.

## Historical Three-layer evidence (50 seeds/cell)

All four source-calibration/target-anchor variants now have300 matched results
at p100/200, C1–C3 and50 seeds/cell. The final report is scratch
`implementation/r8/calibration_ablation_final_review_v1/`. The six p10 execution
checks and600 new high-dimensional results are complete; other factorial cells
reuse checked frozen results. Data and ordinary-target references match.

Source calibration's finite-sample MSE effect is mixed. Omitting it raises
joint-TATE MSE by1.20e-5 (MCSE5.12e-6) at C2/p100, but changes it by-2.25e-6
(MCSE7.24e-6) at C2/p200. Replacing the Hou anchor with ordinary target fitting
raises MSE in all six cells, with substantial Monte Carlo uncertainty. The
experiment is complete; a universal calibration advantage is not established.

## Verified

- The default bounded DGP supplies `W_outcome=X` and `Z_site=X`. The active
  score-derivative solver checks equality of the working bases. C2 changes the
  true outcome feature, while C3 changes the true assignment feature; both
  fitted nuisances still use the same X. Truth misspecification is not a
  different working-basis choice.
- C2 target propensity is explicitly expit(X theta), with
  theta=(.35,-.25,.125,-.0625,0,...). C3 deliberately makes assignment nonlinear
  while its outcome working model is correct. C4 is the explicit both-wrong
  boundary. Thus C2's target anchor is not doubly misspecified by this DGP.
- Calibrated target-only is stored as `target_anchor_ate`. The p100/rho0
  C1–C3 cells each have200 repeats, with coverages96%,96%,95%. Ordinary target
  coverages are94.5%,94%,95%; their RMSEs are very similar. These differences
  do not establish a universal calibration advantage.
- Saved Overleaf revision8be3d007 has the main-text lemma
  `lem:main_calibration` at main.tex:346, its explanatory proof paragraph,
  and supplement `lem:pairwise_orthogonality` at supplemental.tex:621.
  This verifies the lemma's presence, not complete alignment of every theorem
  with the current three-level, fixed-cutoff implementation.

## Remaining manuscript alignment

The same saved manuscript's simulation section still describes the older
reference DGP and its induced target propensity, explicitly noting possible
target double misspecification in Scenario2. It does not describe the current
bounded default or current simulation results. The user's later online edits
have not been freshly synchronized in this check.

At the initial audit, there was no completed full “with/without final score
calibration” comparison. The factorial report above now supplies that comparison.
`target_lasso` changes the target anchor, `initial_validation` changes the
inner source-validation program, and `two_round` changes initialization and
communication. None is the requested complete source-calibration ablation.

## Controlled comparison and historical rollout

The completed controlled two-by-two comparison is:

| Source nuisance program | Target anchor | Status |
|---|---|---|
| Final calibrated fit | Hou calibrated | Complete; current full method |
| Final calibrated fit | Ordinary lasso | Complete; target_lasso ablation |
| Standard outcome fit plus initial merged-weight fit | Hou calibrated | Complete; 300 matched results |
| Standard outcome fit plus initial merged-weight fit | Ordinary lasso | Complete; 300 matched results |

The following paragraphs retain the implementation and interim rollout history.
Their pending-job statements are historical and are superseded by the completed
report above.

Implementation update: the two missing programs now exist behind
`calibration_control$source_nuisance_method="standard"`, with complete inner
validation and either target anchor. The reviewed changes passed 998 R
assertions and 29 workflow tests, including source training exclusions and
same-seed reuse. A C2 default-regression fixture reproduced the previous
estimate, SE, weights, source estimates and outer coefficients exactly.
The frozen release and experiment preparation are under scratch
`implementation/r8/`. This is implementation evidence; the two factorial
cells remain incomplete until their paired Monte Carlo runs are available.
See [the ablation design](calibration_ablation_design.md).

The six p10 production-grid pilots subsequently completed. Source candidate
estimates match across target-anchor choices, data hashes match, ordinary
target estimates/SEs/truths match existing baselines within 1e-12, and saved
inner/outer training exclusions pass. The p100 controllers are submitted;
p200 is gated on p100 first-seed execution checks. These pilots do not provide
a Monte Carlo coverage comparison.

The p100 ordinary-target branch then completed all150 repeats. With
`joint_tate`, standard/calibrated source coverage was90%/90%,90%/90%,96%/96%
in C1–C3. Standard minus calibrated paired MSE differences were
4.996e-6,1.188e-5,-8.880e-6, with MCSE6.484e-6,5.060e-6,4.409e-6.
Thus source calibration's finite-sample TATE effect is small and mixed in
this branch. Do not claim uniform improvement or change the default DGP in
response. The calibrated-target factorial cells are still computing.
The paired report is under scratch
`implementation/r8/calibration_ablation_review_interim_v2/`.

Hold data, outer/inner fold indices, aggregation rule, cutoff2, nuisance grid,
and numerical tolerances fixed. Each inner validation fit must evaluate the
complete selected nuisance program on its allowed training subset. Standard
source outcome fitting should use source outcomes on those same allowed rows;
simply retaining a target initial OR would confound calibration with removal
of source outcome fitting. Initial merged-weight estimation already uses a
balancing loss, so name the comparison “without final score calibration,”
rather than claiming that all balancing has been removed.

Use a clear source-nuisance-program argument and preserve compatibility of
existing controls; do not overload `source_validation_method='initial'`.
The existing validation name may need a canonical “complete program” value
when both calibrated and standard nuisance programs are supported.

Begin with transportable C1–C3 cells to isolate calibration. Retain calibrated
and ordinary target-only comparisons. Report both point-estimation metrics
and inference scope; the standard uncalibrated score need not have the same
orthogonality guarantee under single-model misspecification. Add the optional
effect-heterogeneity panel only as a declared sensitivity experiment, keeping
the default and all frozen campaigns intact.

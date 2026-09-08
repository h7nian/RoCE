# Target-only TATE: recovered states and actual nuisance-system probe

No primary method, original output, working basis, clipping bound or nuisance
tuning rule is changed. This prepares an estimator-level analytic repair;
it is not bootstrap or SE inflation.

## State recovery

The earlier 75-call replay retained only warning-producing model objects.
A storage-only variant under
`results/direct_tate_mc500_b5000/target_nuisance_state_capture_v19/` repeats the
same CV calls and retains every model. Fitting arguments, seeds, prediction
checks and frozen v19 package requirements are unchanged.

Job 18432073 completed 0:0 in 1:45. All 75 calls have zero saved-prediction
error and zero warnings. There are 25 shared propensity states and 50
arm-specific outcome states, covering outer and inner target-only fits.
All output hashes pass. These are replayed states with prediction parity,
not bytewise verification of original model objects that were never saved.
Original CV partition-vector identity remains explicitly unverified, as in
the preceding replay. No new independent MC observations were created.

The model keyed `PS_1_k1_k2` is shared by BOTH arms: its prediction is P(A=1).
It must not be duplicated as an independent control-arm nuisance.
Captured states are under `captured_states/`.

## Actual 603-dimensional system

`probe_target_nuisance_system.R` checks seed10013, outer fold1: one
201-coordinate propensity model and two 201-coordinate outcome models.
Moments are the logistic score for A and arm-specific outcome scores,
normalized by the full target training size. The Jacobian is block diagonal;
shared observations still induce cross-covariance between moment contributions.

The point score retains propensity clipping [.01,.99] and outcome clipping
[.001,.999]. Nuisance loss moments use raw logistic probabilities, so their
derivatives differ from clipped point-score derivatives. The propensity
gradient contains BOTH TATE arms with their correct signs.

Job 18432613 completed 0:0 in 15 seconds; output `target_system_probe_v1/`:

- Saved-prediction error: zero; target fold TATE error: 2.776e-17.
- Twelve joint directional checks: maximum Jacobian error 2.414e-11 and
  score-gradient error 9.175e-12.
- Native-fit KKT checks include glmnet feature standardization and the
  outcome-arm/full-target normalization. Maximum errors: PS 2.763e-5,
  treated outcome 1.277e-7, control outcome 2.398e-5.
- Training clipping counts are zero; held-out propensity clipping occurs
  once. This is distinct from source inference-tilt clipping.

## Coefficient probes, not a selected repair

All three fractions solve their projection equations. Held-out records are
never used to choose a fraction:

| Fraction | Coefficient L1 | Original fold TATE | Adjusted fold TATE | Original score SD | Adjusted score SD |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 0.1677125 | 0.1677125 | 0.9224372 | 0.9224372 |
| 0.5 | 0.82982 | 0.1677125 | 0.1712448 | 0.9224372 | 0.9594072 |
| 0.1 | 29.27384 | 0.1677125 | 0.3244776 | 0.9224372 | 1.3053915 |

The weak-penalty target shift is large. Neither it nor the middle fraction
is adopted from this fold's behavior. These are single-fold diagnostics, not
bias/coverage estimates or a full revised estimator. The retained outer/inner
states permit training-only projection validation without refitting those
nuisances. A rate-valid policy, no outer-evaluation leakage through nuisance
models, source-graph integration and final common weights remain required.

Hashes:

- Capture: `87a5a473b7b234bd83dc91a7c9cf263c44370b9f940bac6a6a7cca28f1972a73`.
- Target probe: `5914cb2fa038e9f72e0446c19452088ea02fa9551fd53fbc49d6da4960bada44`.
- R launcher: `06fd33380baae937240e8b9cbe541bb6d379aa68513b7b8a13c822e238fd3bec`.

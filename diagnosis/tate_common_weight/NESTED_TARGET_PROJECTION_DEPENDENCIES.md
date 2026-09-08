# Inner target projection: exclusion requirements

The five-fold outer projection check does not supply held-out coefficients
for the inner target scores used to learn common TATE weights. Inspection of
`R/cross_fitting_algorithms.R` confirms that each original inner target model
excludes both the outer fold and its inner evaluation fold.

The executable `plan_nested_target_projection.R` audited all 20 ordered
outer/evaluation pairs against the saved five-fold selection states. Every
inner evaluation observation contributes to one validation loss and three
coefficient-training subsets in that outer selection. Reusing its selected
scale would therefore not preserve the same inner holdout property. This
does not prove an existing RoCE bug: outer evaluation is still excluded, and
other theoretically justified weight-learning designs may allow such reuse.
It does rule out describing the naive reuse as fully inner-held-out tuning.

To apply the same validation-risk policy within an inner training complement,
its validation nuisance model must exclude outer fold k, inner evaluation
fold l, and validation fold r. Cached one-/two-fold-excluded fits cannot meet
this three-fold exclusion. The new split manifest contains 60 ordered
(k,l,r) rows but only 10 distinct training subsets. With one shared propensity
and two arm-specific outcome models per subset, a deterministic shared-fit
policy could require 30 models rather than 180. That policy and its seeds
must be frozen before fitting; models from distinct excluded subsets must
never share a cache key.

The smaller training subsets contain only two of the original five folds.
CV class counts, clipping, convergence, projection stability and full score
covariance must consequently be audited again. This is not permission to
discard failures or select a policy by the known TATE.

Evidence is retained in
`results/direct_tate_mc500_b5000/target_nuisance_state_capture_v19/nested_target_dependency_audit_v1/`.
No new models were fitted. No source correction or final analytic inference
is established by this dependency audit.

The subsequent fit policy is fixed in `NESTED_TARGET_FIT_PROTOCOL.md` and
implemented by `fit_nested_target_states.R`. Job 18438592 was submitted for
the 30 shared fits, output `nested_target_states_v1/` under the same root.
Submission is not evidence that fits passed; inspect Slurm status and the
committed diagnostics before consuming the states.

Job 18438592 subsequently completed with exit code 0 in 25 seconds. All 30
fits were eligible under the fixed numerical rules, with zero captured
warnings. `audit_nested_target_states.R` independently checked the exact
training/model IDs for all ten triples, all 30 seeds and balanced CV label
counts, CV training class support, zero glmnet error codes, and coefficient
prediction identity on all 1,000 target observations (maximum error zero).
All five fit-bundle payload checksums passed. The audit bundle is
`nested_target_states_audit_v1/` under the same root.

Two-fold training subsets have 398--402 observations, including only 158--160
treated observations. Some selected treated outcome fits are intercept-only.
For excluded triple 1:2:4, the minimum raw propensity prediction over the full
target dataset is 0.004364264, below the estimator's 0.01 floor. This is not
itself a failed fit: projection validation must use the actual clipped score
and corresponding active clipping derivatives, and distinguish training from
held-out clipping. No positivity or nuisance-rate guarantee follows from the
successful numerical audit.

The shared projection selector now accepts an explicit three- or four-fold
validation set, while preserving the candidate grid, weighted risk and
all-fold failure policy. Its 34 assertions pass. Inner projection selection
has not yet been run with these new states.

The actual inner-score driver is now `run_nested_target_projection.R`, with
the pre-run policy in `NESTED_TARGET_PROJECTION_PROTOCOL.md`. Job 18438695
was submitted for all 20 ordered pairs and 240 candidate/validation attempts,
using cached nuisance states without refits. Its intended output is
`nested_target_projection_v1/`. Completion and downstream suitability must be
determined from the actual job and payload, not this submission record.

Job 18438695 completed with exit code 0 in 37 seconds. All 240 candidate
attempts and 20 selected final coefficient fits succeeded. Ten pairs selected
scale 1 and ten selected null; no other scale was selected. All original
inner prediction reconstruction errors were zero. All four payload checksums
passed. Independent recomputation from the saved per-pair losses reproduced
all 20 selections and verified exclusion of each evaluation ID from its
coefficient-training and validation sets.

The largest selected changes warrant further diagnosis, not immediate
adoption: pair (outer 2, evaluation 4) changes 0.2717830 to 0.3304165 with
score SD 1.2067057 to 1.3724962; pair (outer 4, evaluation 2) changes
0.3063974 to 0.3569216 with SD 1.1766174 to 1.3657274. The corresponding
coefficient L1 norms are 1.6359860 and 1.5650122. Full per-observation scores,
candidate fits and clipping records remain available for decomposition.
No common weights were recomputed; neither these target-only shifts nor the
successful numerical tests establish a valid full-estimator variance.

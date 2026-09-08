# Saved projection: outer-fold score stability

This checks the same three numerical probe coefficients on outer fold1 of
seed10013/rho0/source s1/treated arm. It does not select a penalty, fit a new
model, change a published result, or assess TATE coverage. The original
observations, outcomes and coefficients are retained. No bootstrap or SE
inflation is used.

`audit_saved_projection_holdout.R` verifies the saved matrix and seed artifact
hashes, reloads the pinned v19 package, reconstructs the exact fitting state,
and confirms that every evaluation ID is disjoint from all corresponding
fitting, calibration and projection rows. This is a record-exclusion check,
not a proof that treatment-stratified folds are probabilistically independent.

The held-out moment function includes all ten nuisance blocks. Target final
tilting moments average the four initial outcome derivatives equally; source
final moments use the original calibration-fold size fractions. Initial
tilting moments retain their untruncated exponential loss weights, whereas
the original inference score retains its M=5 truncation. Contributions are
summed on each shared observation before inspecting variability; the initial
fits are not treated as independent samples.

## Results

Execution completed successfully and the checksummed output is
`results/direct_tate_mc500_b5000/independent_inference_pilot_v19/p100_projection_holdout_v1/`.
All payload hashes pass. The original held-out source-assisted arm mean
matches the saved native value, 0.7114726. Collapsed per-observation corrections
and direct block-moment means agree within 4.164e-17. Both held-out sites have
201 observations in this particular outer fold.

| Probe fraction | Coefficient L1 norm | Adjusted fold arm mean | Target score SD | Source score SD | Source max absolute centered score |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 (zero correction) | 0 | 0.7114726 | 0.1128433 | 0.8388452 | 4.468478 |
| 0.5 | 1.39324 | 0.6972403 | 0.1599982 | 0.8395661 | 4.473041 |
| 0.1 | 44.46110 | 0.8222468 | 0.8757095 | 3.0217618 | 24.885622 |

The weak-penalty case increases source score SD about 3.6-fold and target
score SD about 7.8-fold in this fold. Its largest source observation accounts
for 33.91% of the source centered sum of squares. These are held-out function
stability warnings despite successful training KKT checks. The corrected
fold mean is a single-arm/source diagnostic, not the final aggregated TATE.

The weak-penalty mean adjustment decomposes into +0.1270757 at target and
-0.0163016 at source. The initial tilting blocks have substantial individual
target/source terms with partial cancellation; the final outcome block alone
contributes +0.0875417 at source. Coefficient L1 norms include 18.4069 for the
final outcome block and 11.6063 for the final tilting block. Thus the observed
growth is not merely a solver stopping error or an ignored initial block.

All three fractions remain reported. The middle fraction is NOT selected
because its fold mean or variability looks more favorable; neither truth nor
coverage was used for parameter selection. This one diagnostic fold cannot
establish a regularization policy, sampling stability, or improved inference.

## Consequence for the structural repair

Do not promote the raw weak-penalty correction. A training-only, rate-justified
projection policy and control of the resulting score functions are still
needed. Initial/final nuisance coupling and held-out feature leverage matter,
not just the norm of the adjoint equation residual. The target-anchor and
inner aggregation graphs must also be handled. The separate C1 common-weight
stability problem remains unresolved. None of the full primary, sensitivity
or RHC completion requirements is replaced by this check.

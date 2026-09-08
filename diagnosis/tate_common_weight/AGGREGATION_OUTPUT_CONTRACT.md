# Interpreting common TATE weights under cross-fitting

The primary estimator forms target-only and source-assisted TATE contrasts,
then jointly learns one common source-weight vector **per outer fold** using
only its training complement. It does not separately weight the two arms.

## Exact finite-sample reconstruction

For outer fold k, write T_k for `fold_target_estimate`, M_jk for
`mu_pred_ts[j]`, D_jk for `delta_ts[j]`, and eta_jk for `fold_weights[k,j]`.
Both M and D are already treated-minus-control contrasts. With target fraction
p_tk = n_tk/n_t and source fraction p_jk = n_jk/n_j, the exact estimate is

```text
sum_k { p_tk * [(1 - sum_j eta_jk) * T_k + sum_j eta_jk * M_jk]
        + sum_j p_jk * eta_jk * D_jk }.
```

This is the mean of the site-scaled observation-level pseudo-values in
`.compute_phase3_all_phi()`. When every site has the same fold fractions,
it reduces to the corresponding weighted average of fold-level aggregated
TATEs. Using the actual fractions avoids assuming perfectly balanced folds.

Outer `intermediates$fold_info[[k]]$source_idx` lists may be unnamed. Their
positions follow `names(weights)`, as checked by the saved-fit validator.
Resolve `j = match(site, names(weights))` and read `source_idx[[j]]`; do not
assume that `source_idx[[site]]` returns IDs. Reject a missing source position
or empty/noninteger/out-of-range IDs before forming a training complement.
Otherwise an accidental NULL lookup can silently count the entire source as
training data. This is an observation-index interpretation rule, not a change
to the estimator or a new set of source weights.

## Reporting fields are not a second aggregation step

`fold_weights` records the weights actually used. The returned `weights`
vector is their arithmetic mean, useful as a source-borrowing summary.
`source_estimates` pools the source-assisted estimates with the appropriate
target/source fold fractions. In general,

```text
estimate != target_only$estimate +
            sum(weights * (source_estimates - target_only$estimate)).
```

Averaging weights and estimates separately loses their foldwise association.
The right-hand side instead describes a different estimator that uses that
constant vector in every fold. No unique set of effective global source
weights is implied by the cross-fitted result. Do not invent such weights
just to make the pooled display agree with the estimate.

The six saved v18 fits for seed 1 give the following pooled-display minus
actual-estimate differences. These are reporting-reconstruction discrepancies,
not Monte Carlo bias and not evidence that the implemented aggregation is wrong.

| rho | Mean-weight pooled display minus actual estimate |
| ---: | ---: |
| 0 | -0.00661795307 |
| 0.5 | -0.01220059409 |
| 1 | -0.00889749598 |
| 1.5 | -0.00324831290 |
| 2 | -0.00324831290 |
| 2.5 | -0.00324831290 |

The exact finite-sample reconstruction above agrees with all six actual
estimates to at most 5.6e-17. `test_full_refit_resampling.R` checks it and the
constant-vector pooled identity against saved artifacts.

## Variance-derivative bookkeeping

Relative to any reference weights eta*_jk, the exact change in the estimator is

```text
sum_k sum_j (eta_jk - eta*_jk) *
  [p_tk * (M_jk - T_k) + p_jk * D_jk].
```

Only when p_jk = p_tk does the bracket simplify to p_tk times the fold-level
source-target TATE discrepancy. This identity is algebraic; whether its terms
are asymptotically negligible requires separate rate and stability arguments.

Two regimes must be distinguished:

- For a nonzero limiting source-target discrepancy, root-N variation of an
  active learned weight can contribute a first-order term proportional to
  that discrepancy. In a regular differentiable regime its influence is of
  the form `b^T IF_eta`, including covariance with the fixed-weight influence.
  This expression is not automatically valid at selection boundaries.
- Even when the limiting discrepancy is zero, a nonconvergent random weight
  can multiply an O_p(N^-1/2) empirical discrepancy and remain first order.
  Zero population discrepancy alone is therefore not enough to discard
  weight uncertainty. Consistency toward fixed reference weights plus the
  compatible-source expansion is a sufficient negligibility route.

The manuscript's oracle route excludes incompatible sources with probability
tending to one and establishes weight consistency on compatible sources under
a rate-indexed penalty sequence. Together with the stated nuisance expansion
and variance regularity conditions, this supplies the missing negligibility
argument. `docs/main.tex` explicitly restricts that theorem to lambda_N tending
to zero while lambda_N sqrt(N) diverges; it does not automatically validate the
implemented fixed cutoff c=1. A cross-fitting training/evaluation split alone
does not prove unconditional validity for the complete adaptive estimator,
because observations also influence weights used in other folds.

The expanded saved-fit contract tests passed 143 assertions on seed 1's six
v18 artifacts. This verifies the finite-sample algebra, not the asymptotic
conditions or finite-sample coverage.

At the next package/documentation freeze, clarify the `weights` reporting
contract in public result documentation. Do not rename fields or modify
frozen package bytes while the current calibration jobs are running.
Also qualify the unconditional-validity wording in the roxygen Phase 2
description of `calculate_crossfit_aggregation()`: independence of the held-out
fold is necessary bookkeeping, not the complete adaptive-weight argument.

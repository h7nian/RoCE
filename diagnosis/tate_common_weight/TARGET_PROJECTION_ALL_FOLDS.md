# Five-fold extension of the target projection development check

The original coefficient policy in `TARGET_PROJECTION_VALIDATION_PROTOCOL.md`
is unchanged. The driver now accepts an optional outer fold (1 through 5),
defaulting to 1 for compatibility. For outer fold k, validation uses precisely
the other four folds; each coefficient-training subset excludes both k and
its current validation fold. Model keys explicitly include k. The selector
receives the validation fold set instead of assuming folds 2 through 5.

Job array 18437974 submits these five checks, with at most two concurrent
single-CPU tasks, using the frozen v19 package and cached seed 10013 fits.
Outputs are separate `target_projection_fold_<k>_v2` directories under
`results/direct_tate_mc500_b5000/target_nuisance_state_capture_v19/`.
The original single-fold bundle and its source snapshots remain untouched.
The revised validation suite passed 29 assertions before submission.

Every task must pass prediction parity with its saved outer and inner fits,
record coefficient-fit failures without substituting another candidate, and
commit its payload checksums. This is not a new independent replication,
coverage gate, source correction, or full fitted-estimator variance estimate.
Any later aggregation must restore observation IDs and use actual fold sizes,
not an unweighted mean of fold estimates.

## Completed check

All five array tasks completed with exit code 0 (9--17 seconds each).
`audit_target_projection_folds.R` verified all six payloads per fold, identical
fitting-code hashes across folds, reproduced candidate risks and selections,
and checked all outer/inner exclusion sets against the complete 1,000 target
observation IDs. It restored the original observation order and verified that
the pooled means equal the actual-fold-size-weighted fold means.

Folds 1, 2, 3 and 5 selected penalty scale 1; fold 4 selected the null
projection. No fit failed and no fold was discarded. The target-only pooled
estimate changed from 0.1670652 to 0.1745530. Fold 2 had the largest change,
0.2692002 to 0.2971062, and score SD increased from 1.0701292 to 1.1704045.
These are descriptive results from one development replication, not evidence
of improved RMSE or coverage. In particular, larger score SD is not itself a
valid correction to the full fitted-estimator standard error.

The audit output is `target_projection_five_fold_audit_v1/` under the same
root; all four output payload checksums passed. It deliberately reports no
confidence interval and sets `inference_validated=FALSE`. The next integration
must address the inner target scores used in weight learning and the source
nuisance graph together; replacing only the outer target score would leave
the full estimator and weight objective misaligned.

# Source coefficient selection: staged development protocol

Before examining all-fold candidate rankings, fix the following candidate
procedure at seed 10013, outer fold 1, source s1. This is informed by prior
development checks on that seed and is not independent statistical validation.

Stage 1 compares a null final projection and final penalty scales {0.5,1,2}.
For each validation fold, final coefficients are fitted using only its audited
training-complement nuisance state. Sum the two final-block conditional risks
over both arms; these final equations have no dependence on initial projection
coefficients, so this is a valid block-diagonal final-stage loss. Average folds
with their actual target-plus-source validation sample counts. Select minimum
risk, with numerical ties preferring null and then the larger scale. A candidate
is ineligible if any required arm/fold fit fails; never average only successes.

Stage 2 compares null initial correction and initial scales {0.5,1,2},
conditional on the selected stage-1 scale. Enumerate all stage-1/stage-2
branches before selection. Within a branch the downstream final coefficients
are fixed training-fitted vectors; sum initial-block conditional risks over
all initial blocks and both arms, then weight validation folds as above.
Use the same complete-fold failure and tie rules. Do not compare upstream
risks across different downstream branches. If stage 1 selects null, the
coupling right-hand sides vanish and the full null projection is retained.

The final coefficients used within a validation fold are always those fitted
on its complement, not coefficients refitted after looking at that fold.
The selected scale itself depends on validation losses; the two-stage procedure
is therefore adaptive. Finite candidate enumeration preserves the data-flow
audit, but it does not by itself prove nuisance rates or valid Wald inference.

After both choices, refit on the complete outer-training graph using its
original cached nuisance state, with the selected final and initial scales.
Do not inspect the outer evaluation records until those choices are fixed.
No TATE truth, coverage, interval width, or post hoc clipping change enters
selection. Repeatable sampling and full target/source/weight covariance
validation remain required before adopting this method.

The constructor retains its old single-scale default. Passing an explicit
`initial_penalty_scale` separates the stages; NULL denotes a no-initial-
correction candidate and must not be labeled a solved unpenalized adjoint
equation. Retain all such candidate labels and residuals.

## Selector implementation check

`source_stagewise_selection.R` requires the complete 128-row table (four
final candidates, four initial candidates, two arms, four folds). Final risks
must agree across initial branches, and validation counts must agree within
fold. It sums arms before applying the shared complete-fold selector, keeps
final/initial failures separate, and conditions stage 2 on the selected final
branch. All 13 unit assertions pass, including branch isolation, row-order
invariance, missing rows, inconsistent risks, a failed initial candidate, and
one failed final arm. This is a selector software check, not an executed
all-fold source selection or statistical validation.

The executable `run_source_stagewise_selection.R` now enumerates the full
candidate table, retains failures by stage, invokes the selector and refits
the selected coefficient rules on the original outer-training state before
evaluating outer fold 1. Job 18443364 depends on successful completion of
all audits in array 18442608. Intended output is `source_stagewise_selection_v1/`
under the independent-pilot root. This job performs no new nuisance fits and
does not recompute common TATE weights or claim validated intervals.

The driver now accepts an explicit C1/C2/C3 configuration (default C1).
C2/C3 selection reads matching basis-specific validation audit bundles and
matching complete-outer-training audit states; it never substitutes C1
coefficients. Basis dimensions and exact full-state reevaluation are checked.
C1 regression job 18444629 writes `source_stagewise_selection_v2/`; inspect
its scientific outputs against v1 before claiming backward compatibility.

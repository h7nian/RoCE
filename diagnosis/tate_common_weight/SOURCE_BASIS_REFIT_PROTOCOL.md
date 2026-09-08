# Paired C2/C3 working-basis refits

Use seed 10013's fixed FACE observations at rho 0, source s1, outer fold 1
and validation fold 2. `source_working_basis.R` follows the current data
generator's working-basis definitions: C2 uses raw X for outcome regression
and retains [X,X^2] for site calibration; C3 retains the original centered
quadratic outcome basis and uses raw X for site calibration.

The transformer verifies the input quadratic blocks. On all three saved
sites, checks confirmed that X, X_dagger, A, Y, sample sizes and true-model
fields are unchanged. C1 is an exact identity operation. Fits use the same
seeds and nuisance settings as the C1 first validation split, differing only
in working basis. This is a paired diagnostic, not newly sampled C2/C3 Monte
Carlo replications or completed configuration-level coverage evidence.

Fit array 18443433 covers C2 and C3; its output roots are
`source_basis_C2_split_1_2_v1/` and `source_basis_C3_split_1_2_v1/` under
the independent-pilot directory. Equation audits use the saved working_config
to reconstruct the same basis. Older C1 bundles without this field default
to C1. Native fits, failures and warnings must be retained. Numerical success
does not by itself resolve the known calibration tangent-span issue or prove
double robustness; full source projection and inference require further checks.

The dedicated `test_source_working_basis.R` suite passes 15 assertions:
C1 identity, correct C2/C3 changed field, unchanged observations and truth,
no mutation of the input, invalid configuration rejection, and rejection of
already reduced or malformed quadratic inputs. Audit array 18443436 depends
on fit array 18443433. No runtime or numerical completion claim is implied by
these transformation tests.

Both paired fits completed successfully: C2 in 6:11 and C3 in 13:02, with
zero captured warnings in both arms. Their audits (18443436_2 and _3)
completed in 12 and 14 seconds. All fit and audit payload checksums passed.
Each arm's unequal-basis system has 1,208 parameters. Across all four arm/
configuration systems, maximum Jacobian error is 4.271e-10, score-gradient
error 7.407e-11, final KKT error 3.618e-7 and held-out point error zero.
These establish actual-fit equation consistency, not double robustness or
coverage. Coupled projection behavior under these fits remains to be checked.

After the first-split candidate check, fit array 18443702 was submitted for
the six remaining configuration/validation pairs: tasks 1--3 are C2 folds
3--5, tasks 4--6 are C3 folds 3--5. At most two single-CPU tasks run at once.
All use the existing candidate-independent nuisance settings and the same
fold-specific seeds as the C1 counterparts. No failed or unfavorable split
may be omitted from subsequent selection. These remain fixed-observation
paired diagnostics, not additional independent replications.

Audit array 18443708 uses the identical task-to-configuration/fold mapping
and depends on successful completion of 18443702. The full outer-training
C2/C3 nuisance states are not present in the original C1 seed bundle: they
must be fitted and audited separately before any selected outer-score
evaluation. Do not combine a C2/C3 validation choice with C1 full-training
coefficients merely because the underlying observations are shared.

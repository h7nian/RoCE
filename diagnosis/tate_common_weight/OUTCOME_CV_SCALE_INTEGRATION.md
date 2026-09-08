# Tested normalization correction integrated into workspace source

Before editing, workspace `src/outcome_model.cpp` was byte-identical to the
frozen v19 source copy, although it contains existing user/worktree changes
relative to git. Those changes were preserved. Only the two tested CV training
scale expressions and their comments were patched in the workspace.

`tests/testthat/test-outcome-cv-training-scale.R` now contains the same 16
assertion regression used in the isolated candidate. The root copy passes
against the verified installed candidate. Workspace DESCRIPTION, NAMESPACE,
all R files, C++ sources/headers and Makevars were compared byte-for-byte to
the software-tested candidate and all match. The preceding candidate evidence
includes `Status: OK`, 2,293 repository assertions, 16 targeted assertions,
and seven separately enabled related audit assertions.

The frozen v19 installed/source directories and all archived experiments are
unchanged. No historical result is relabeled as having the corrected CV scale.
This workspace correction does not adopt either experimental source projection
or shared-final-basis estimator, and does not resolve the full TATE inference
problem. Final production freeze and experiments still require coherent
source/package/workflow hashes and statistical validation of the final method.

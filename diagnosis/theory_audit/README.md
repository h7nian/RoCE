# Independent theory audit

This directory implements the first, theory-audit stage of the agreed repair
plan. Production code and the manuscript are unchanged. The user selected a
theory review before deciding changes to the core method; the findings require
that decision before estimator redesign or new performance production runs.

Report, claim ledger, function-review ledger and all outputs:

[Theory report](/scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/theory_report.md)

The principal report is `theory_report.md` in that directory. Its twelve
findings distinguish counterexamples, proof gaps and implementation scope.
Population fixtures are not estimates of FACE simulation coverage.

The user designated `/scratch.global/zhan9381/FACE-HD/` as the canonical result
root. Earlier audit trees were moved there with atomic renames; old paths are
compatibility symlinks for recorded provenance. See `diagnosis/storage_migration.json`
under that root. New files and temporary outputs use the canonical root.

`hou_comparison.md` records the subsequent paper review. The accompanying
`hou_calibration_checks.py` checks three model-correctness cases for both arms,
including distinct initial/final limits. It verifies population equations, not
the complete SMMAL estimator or its high-dimensional theorem.

## Reproduce

Run from the project root. Each command requires a **new** output directory.
Python needs NumPy and SciPy. R inventory uses base R and does not load project
code; the package fixture uses the explicitly supplied frozen library.

```bash
python3 -B diagnosis/theory_audit/build_audit_index.py --output /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/source_index_new
OPENBLAS_NUM_THREADS=1 python3 -B diagnosis/theory_audit/population_checks.py --output /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/population_checks_new
OPENBLAS_NUM_THREADS=1 python3 -B diagnosis/theory_audit/hou_calibration_checks.py --output /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/hou_calibration_checks_new

module load R/4.2.2-gcc-8.2.0-vp7tyde
TMPDIR=/scratch.global/zhan9381/FACE-HD/tmp Rscript --vanilla diagnosis/theory_audit/inventory_r_functions.R /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/source_index_new/source_manifest.csv /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/r_function_index_new
R_LIBS_USER=/users/0/zhan9381/Rlibs ROCE_PROJECT_LIB=/users/0/zhan9381/FACE-HD/results/direct_tate_mc500_b5000/Rlib_production_20260917_v4 TMPDIR=/scratch.global/zhan9381/FACE-HD/tmp OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla diagnosis/theory_audit/check_package_population.R /scratch.global/zhan9381/FACE-HD/diagnosis/repair_audit/package_population_new
```

The Python checks independently verify ten population identities or
counterexamples. The package fixture separately checks the actual nuisance
primitives with exactly 1000 observations per site at fixed zero penalty.
It does not run the complete cross-fitting/aggregation pipeline.

R functions include nested and anonymous definitions. Parameter body mentions
and assignment names are navigation aids, not tests of parameter semantics or
evidence for deleting code. Python definitions use the AST; C++ and shell
definitions are explicitly labeled lexical candidates. Slurm `.cmd` launchers,
build configuration and existing generated Rcpp interfaces are included.
Historical results, logs, the Overleaf copy and other git-ignored files are
excluded from source discovery.

`claim_review.csv` and `function_review.csv` are reviewed annotations, separate
from generated discovery indexes. Pending semantic reviews remain pending;
neither a source-text match nor passing syntax checks certifies correctness.

# Reproduce the current RHC analysis

For All of Us or another dataset, use [the general real-data guide](REAL_DATA.md).
RHC's public teaching CSV is included in the installed package.

After the installation instructions in that guide:

```bash
export ROCE_CORES=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript scripts/real_data/run.R scripts/real_data/rhc.R /path/to/output/rhc preflight
Rscript scripts/real_data/run.R scripts/real_data/rhc.R /path/to/output/rhc all
```

On MSI an example output root is
`/scratch.global/zhan9381/FACE-HD/real_data/rhc/collaborator_reproduction/`.
Use your own writable location elsewhere. RHC exact reproduction expects the
Linux/WSL `C.UTF-8` collation used by the frozen analysis.

The explicit profile keeps the Private-insurance target (1698 patients), three
sources (Medicare 1458, Private & Medicare 1236, Medicaid 647), 5039 patients in
total, 61 working features, death within 30 days, historical ten-fold assignment,
Two-layer cross-fitting, one-round communication, joint TATE, cutoff 2, original
`min` nuisance CV with 100 penalties, source truncation radius 3 and target radius 5.
The two excluded insurance groups remain excluded. No target rotation or
source-1se alternative is part of this profile.

Imputation, centering and scaling exclude all sites' current outer evaluation
fold. Initial/calibration CV shares that outer-training transformation. SS/IVW
use their original site-specific folds and exclude the site's validation rows
from preprocessing; Federated-DR/Pooled-DR retain full-sample fitting. Baselines
use paired-arm influence covariance and 5000 bootstrap multipliers.

The frozen analysis to compare with (outcome scale, not percentage points):

| Method | Estimate | SE |
|---|---:|---:|
| Target-only (calibrated) | 0.0240017838 | 0.0246111620 |
| RoCE | 0.0301579860 | 0.0216427172 |
| SS | 0.0599344470 | 0.0156863585 |
| IVW | 0.0583093805 | 0.0159419496 |
| Federated-DR | 0.0438655238 | 0.0193523358 |
| Pooled-DR | 0.0638409135 | 0.0272665132 |

Small compiler/library differences can change numerical tolerances. Compare
input identity, complete configuration, diagnostics and point estimates before
interpreting SE changes. The SE comparison is descriptive; it does not establish
causal bias or coverage on this dataset. Raw fit objects remain available only
if explicitly requested via `ROCE_SAVE_FIT=1`.

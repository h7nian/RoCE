# Reproduce RHC

After [installation](../README.md#install), run the bundled public RHC data:

```bash
export ROCE_CORES=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript scripts/real_data/run.R scripts/real_data/rhc.R /path/to/output/rhc preflight
Rscript scripts/real_data/run.R scripts/real_data/rhc.R /path/to/output/rhc all
```

Exact reproduction uses Linux/WSL `C.UTF-8` collation. For other datasets,
use the [custom-data guide](REAL_DATA.md).

The [profile](../scripts/real_data/rhc.R) uses 30-day death, 61 features and
5039 patients: Private target (1698), Medicare (1458), Private & Medicare (1236),
and Medicaid (647). Settings are ten folds, two-layer/one-round, joint TATE,
cutoff 2, `min` CV with 100 penalties, source radius 3 and target radius 5.

Preprocessing excludes the current outer evaluation folds; inner CV shares the
outer-training transform. SS/IVW retain their site-specific folds and training
transforms. Federated-DR/Pooled-DR retain full-sample fitting. Baselines use
paired-arm covariance and 5000 bootstrap multipliers.

Reference results on the outcome scale:

| Method | Estimate | SE |
|---|---:|---:|
| Target-only (calibrated) | 0.0240017838 | 0.0246111620 |
| RoCE | 0.0301579860 | 0.0216427172 |
| SS | 0.0599344470 | 0.0156863585 |
| IVW | 0.0583093805 | 0.0159419496 |
| Federated-DR | 0.0438655238 | 0.0193523358 |
| Pooled-DR | 0.0638409135 | 0.0272665132 |

Compare configuration and diagnostics when reproducing these numbers. Real-data
SE comparisons do not establish causal bias or coverage.

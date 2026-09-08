# `diagnosis/c2/` - RoCE C2 coverage diagnosis

Measurement-only C2 experiments submitted through MSI/Slurm.

Use:

```bash
bash diagnosis/c2/c2_diagnose.sh --dry-run
bash diagnosis/c2/c2_diagnose.sh --tasks 1-4 --n-sims 100
```

Task groups:

- `1-12`: exact C2 settings where the main result table showed under-coverage,
  paired as superpopulation/sample estimands at the default `K_f = 10`.
- `13-14`: larger-n K3/p10 C2 checks at the default `K_f = 10`.

Outputs are written to `diagnosis/c2/results/`; logs are written to
`diagnosis/c2/log/`; checkpoints are written to `diagnosis/c2/checkpoints/`.

# Fixed K=4 and K=6 sensitivity experiments

Authorized by the user on 2026-09-21 for the current paper. These are fixed-K
sensitivity studies; they do not establish the next paper's growing-K theory.

| Panel | Dimensions | Scenarios | Outcome deviation | Deviated sources | Repeats/cell |
|---|---|---|---|---|---:|
| Add sources with one incompatible source | 100,200 | C1,C2,C3 | 0,.5,1,2 | 1, with K=4 or6 | 50 |
| Hold the incompatible fraction at one half | 100 | C1,C2,C3 | 1,2 | 2 of4 or3 of6 | 50 |

The first panel has2400 method repeats. The second has600. Each has matching
baseline repeats, giving3000 method and3000 baseline jobs in total. Rho0 does
not change outcomes regardless of the declared number of deviated sources;
it is not redundantly rerun in the fraction-control panel.

Each site retains1000 observations: total sample sizes are5000 for K=4 and7000
for K=6. The DGP definition, sparsity4, ten folds, three-level main method,
one-round protocol, cutoff2, nuisance grid100 and solver tolerance1e-10 remain
unchanged. The stronger-source-count panel is paired within each K by seed.
Changing K changes the generated dataset and the DGP's j/K source-tilt
positions; the study must not be described as adding patients to an otherwise
identical realized dataset. Population TATE remains fixed for each scenario.

Holding only one source invalid reduces the invalid fraction from1/2 at K=2
to1/4 or1/6. The half-invalid panel is a representative p100 high-dimensional
control for that confounding. The p200 fraction-control extension can be
assessed after these diagnostics rather than expanding every interaction at
once. C4 and mechanism ablations continue in their separately authorized K=2
campaigns.

Implementation reuses the existing workers. `--sources` now explicitly
supports6. `--n-deviated-sites` controls the number of leading source outcome
models shifted by rho; its default remains1. Both method and baseline workers
read the same frozen configuration and report the count in their results.
The data-identity audit distinguishes different active deviation counts while
recognizing that rho0 produces no active source deviation.

All outputs are under
`/scratch.global/zhan9381/FACE-HD/implementation/r7/source_count_validation_v1/`.
Source code and readable design notes remain in the current repository.
Checkpointing, immutable input hashes, exact seed ownership and reviewed retry
rules are inherited from the established workflow. No currently running
campaign or checked estimator library is edited.

The new worker-cap sum is432 (288 method,144 baseline). Method repeats request
2 CPUs and8–12GB depending on the profile, with checkpointed two-hour attempts.
Baseline repeats request2 CPUs and4GB. One repeat remains one Slurm job.
First-seed gates check execution and integrity, not desirable statistical
performance. Results include coverage uncertainty, bias/RMSE, SE calibration,
arm-specific source weights and negative-transfer diagnostics.

## Weak-deviation fraction control

The follow-up `source_count_weak_fraction_v1` adds rho=0.5 to the half-invalid
panel at p=100: K=4 with two deviated sources, and K=6 with three. Each uses
C1/C2/C3 and 50 repeats per cell, totaling 300 method and 300 matched baseline
repeats. This separates increasing K from reducing the fraction of incompatible
sources at the weak-deviation boundary. Fixed rho=0.5 is not an asymptotic
local-alternative sequence.

The follow-up reuses the checked scientific configuration and seed rules.
Its separate frozen workflows, results and submission records are under
`/scratch.global/zhan9381/FACE-HD/implementation/r7/source_count_weak_fraction_v1/`.
The additional worker cap is 144 (96 method and 48 baseline). Each method
campaign requires the corresponding earlier half-invalid execution checks;
its matched baselines wait for its own first-seed checks. No coverage or RMSE
threshold determines release or retention of a repeat.

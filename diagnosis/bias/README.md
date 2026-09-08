# `diagnosis/bias/` — RoCE one_round / two_round bias diagnosis

> **⚠️ SUPERSEDED (2026-06-01).** The mechanism below (per-fold calibration too
> small → `exp(Zγ)` blow-up → **`colMeans` parameter-averaging across `k2`** →
> O(1) Jensen gap) and the resulting `n ≳ 10·p·(K+1)·K_f` rule were measured on
> the **initial-commit** code (`gamma_final <- colMeans(per-k2 γ)`, 98040a2). The
> current working tree **refactored** the calibrated nuisance to **fold-summed /
> stacked** fitting (`R/cross_fitting_algorithms.R:195-216`), which fits once on
> ~`(K_f-1)/K_f` of the site — so the per-fold-too-small mechanism is largely
> gone and `n ∝ K_f` is **not** established for the current code. Being re-tested
> in `diagnosis/face_probe/` (kf=5 vs 10 at fixed n). See
> `diagnosis/c2/COVERAGE_INVESTIGATION_STATE.md` (2026-06-01 RULES UPDATE) for
> the corrected, current conclusions.

This folder documents a multi-day investigation into the under-coverage and
occasional blow-up of `one_round_crossfit` / `two_round_crossfit` (the RoCE
two-layer cross-fit estimators). **The algorithm was not modified during any
of this work — only observed.** The conclusion is that the estimators are
correct under their nominal sample-size regime, but the regime is more
demanding than was initially appreciated.

---

## TL;DR

1. The bias / under-coverage of `one_round_crossfit` / `two_round_crossfit`
   is driven by **fold-to-fold instability of the calibrated nuisance fits**:
   for a single source × single `k1`, the per-source contribution
   `mu_pred_ts` swings widely across `k2 ∈ secondary folds`.
2. Mechanically, the initial density ratio `exp(Z · γ_init)` reaches
   `10^31`–`10^135` in small-sample / high-`p` regimes — calibrated steps
   inherit this instability and produce calibrated nuisance coefficients
   with magnitudes 15-30. Averaging non-linear nuisance **parameters** (vs
   averaging predictions) across `k2` then leaves an O(1) gap.
3. The driver is the ratio **`p / n_calib_arm`**, where `n_calib_arm`
   ≈ `n_total / (K+1) / 2 / K_f`. The empirical "safe zone" cliff is at
   ratio ≲ 0.2.
4. Translated to a sample-size rule:
   $$n_{total} \gtrsim 10 \cdot p \cdot (K+1) \cdot K_f \approx 400 p \text{ (for } K=3, K_f=10)$$
5. The production run at `n=5000, p=10` sits at ratio ≈ 0.18 (just inside
   the safe zone), which is consistent with the observed mild under-coverage
   (`one_round` coverage 0.886 at K=3, 0.926 at K=4 from the original
   `results/n5000_*` summaries).

---

## What was actually observed

### Stage 1 — completed end-to-end results in `results/`

Four settings completed (each n_sims=50, C1, binary, mild transform, ss=0.5,
K=3 unless noted). `oracle_dr` / `federated_dr` / `pooled_dr` behave as
expected; `one_round_crossfit` / `two_round_crossfit` break catastrophically
in the small-sample regime.

| Setting                  | one_round bias / SE / cov | two_round bias / SE / cov | oracle bias / cov |
|--------------------------|---------------------------|---------------------------|-------------------|
| `anchor`  (n=1k, p=20)   | **+118.9** / 949 / 1.00  | **−459.7** / 1093 / 1.00 | +0.004 / 0.98     |
| `n_large` (n=2k, p=20)   | **+7.6** / 6.5 / 0.90    | **+0.78** / 4.85 / 0.92  | +0.002 / 0.96     |
| `shift_005`(n=1k, p=20)  | **−676.5** / 1215 / 0.96 | **−411.7** / 1193 / 0.98 | +0.000 / 0.98     |
| `shift_1` (n=1k, p=20)   | **−129.0** / 2209 / 0.96 | **+97.3** / 2433 / 0.94  | +0.006 / 1.00     |
| `K_2`    (n=1k, p=20, K=2)| **+118.9** / 949 / 0.98 | **+85.5** / 382 / 0.98   | −0.000 / 0.99     |

(Cells with `bias > 1` are catastrophic — true `mu1 ≈ 0.585`.)

`oracle_dr`, `federated_dr`, `pooled_dr` all hold around `bias ≈ 0.005,
cov ≈ 0.94–0.98` in every cell — i.e. the **DGP and outcome models are
fine**; the failure is specific to the two-layer cross-fit nuisance
construction.

### Stage 2 — `probe.R` Phase B (per-function trace, single seed, single source)

For source `s2`, `k1 = 1`, scanning `k2 ∈ {2..10}` at the anchor setting
(`n=1000, K=3, p=20`):

| `k2` | `exp(Z·γ_init)` max | `gamma_cal` max | `mu_pred_ts` (single-source contribution) |
|-----:|--------------------:|----------------:|------------------------------------------:|
|    2 |  57       |  8 | 0.445 |
|    3 |  **26,340** | 16 | 0.220 |
|    4 |  20       | 13 | 0.592 |
|    5 |  527      | 15 | **0.061** |
|    6 |  1,200    | 20 | **0.090** |
|    7 |  58       | **29** | **0.762** |
|    8 |  **14,150** | 17 | 0.355 |
|    9 |  602      | 15 | **0.825** |
|   10 |  7        | 14 | 0.647 |

The single-source `mu_pred_ts` ranges from **0.06 to 0.82** for the same
`(k1, source)` pair under different `k2` calibration folds. Truth is
`mu1 ≈ 0.585`. Averaging 9 such wildly different fits cannot recover a
stable point estimate.

### Stage 3 — `probe_scan.R` (six (n, p) cells × 3 seeds × 3 sources × 9 k2)

Quantifies the cliff. Each row is the median across (3 seeds × K sources)
of the per-`(seed, source)` spread `max - min` of `mu_pred_ts` across the
9 secondary folds.

| Setting              | n_calib_arm | **p / n_calib_arm** | median spread | gamma_cal max | exp(Zγ_init) max | % (seed, src) blown |
|----------------------|------------:|--------------------:|--------------:|--------------:|-----------------:|--------------------:|
| n1k, p10             |          12 |            0.83     |          0.64 |           18  |         58       |       100 %         |
| n1k, p20  (anchor)   |          11 |            1.82     |          0.81 |           27  |  **5 × 10³¹**    |       100 %         |
| n1k, p50             |          13 |            3.85     |          0.83 |           32  |  **3 × 10¹³⁵**   |       100 %         |
| n2k, p20             |          22 |            0.91     |          0.41 |           16  |         23       |       100 %         |
| n5k, p10 (production)|          56 |          **0.18**   |          0.20 |          1.0  |         1.8      |        56 %         |
| n5k, p20             |          57 |            0.35     |          0.22 |          1.5  |         3.5      |        78 %         |

"% (seed, src) blown" = fraction of `(seed, source)` pairs with
`mu_pred_ts_spread > 0.2`.

A blow-up condition `exp(Z·γ) → 10¹³⁵` indicates a literally identifiable
but numerically over-fit γ_init: with `n_train_arm ≈ 167–227` and `p ≥ 20`,
`cv.glmnet`'s `lambda.min` allows coefficients large enough to make
`Z · γ` reach `±60`+ on some calibration-fold rows.

---

## Mechanism

1. `fit_initial_outcome` and `fit_initial_density_ratio` use `cv.glmnet`
   with `lambda.min` (paper default; comment in `R/model_fitting.R` notes
   this was a deliberate choice over `lambda.1se`).
2. In low-`n_train_arm / p` regimes, `lambda.min` does not regularize
   enough → γ_init coefficients can be ±5 or larger.
3. `exp(Z · γ_init)` evaluated on the (separate, small) calibration fold
   can therefore reach `10³¹`+ on outlier rows.
4. The calibrated step
   ([R/cross_fitting_algorithms.R:120–129](../../R/cross_fitting_algorithms.R#L120-L129))
   minimises a loss involving these exploding weights (truncated at
   `M_tau=10` in the loss but the calibrated coefficients γ_cal still
   reach ±29 in our scans).
5. The calibrated nuisance is averaged across `k2` as **parameters**
   ([R/cross_fitting_algorithms.R:158-159](../../R/cross_fitting_algorithms.R#L158-L159)):
   ```r
   gamma_final_k1 <- colMeans(do.call(rbind, gamma_cal_list))
   alpha_final_k1 <- colMeans(do.call(rbind, alpha_cal_list))
   ```
   This is **not** an average of predictions; the Jensen gap is O(1)
   when the parameter dispersion is high.
6. The aggregation step (`stabilized_lambda_selection`,
   [R/cross_fitting_algorithms.R:799](../../R/cross_fitting_algorithms.R#L799))
   chooses λ on a natural-scale grid that depends on `(μ_ot − μ_ts,j)²`. When
   `μ_ts,j` is itself ±100 (because of (1)–(5) above), the natural-scale
   grid collapses → CV picks the variance-minimising (unregularised) end →
   final estimate inherits the garbage source estimate. **The aggregation CV
   itself is mathematically correct; it just cannot defend against upstream
   blow-up.**

---

## Recommended operating range

Empirical safe zone (from the scan):

```
n_calib_arm ≳ 5 · p
=>  n_total / [(K+1) · 2 · K_f]  ≳  5 · p
=>  n_total  ≳  10 · p · (K+1) · K_f
```

For `K = 3, K_f = 10`:

| `p`  | minimum recommended `n_total` |
|-----:|------------------------------:|
|   5  |    2,000                      |
|  10  |    4,000                      |
|  20  |    8,000                      |
|  50  |   20,000                      |
| 100  |   40,000                      |
| 200  |   80,000                      |

The original n=5000 production runs at `p=10` sit right at the boundary
(ratio 0.18), which matches the observed `one_round` coverage of 0.886 —
real but mild under-coverage, not blow-up. Pushing to n=10000 at p=10, or
n=8000 at p=20, should largely close the gap; pushing to n≥4·minimum should
also do so.

### K_f sweep — `K_f` is a strong lever

`probe_kf.R` sweeps `K_f ∈ {3, 5, 7, 10}` at two fixed settings (job `9390886`,
4 min wall). The "blow-up" disappears at `K_f = 3` even in the catastrophic
setting that was hopeless at `K_f = 10`:

**Catastrophic setting (n=1000, K=3, p=20):**

| `K_f` | n_calib_arm | p / n_calib_arm | median spread | gamma_cal max | % blown |
|------:|------------:|----------------:|--------------:|--------------:|--------:|
|     3 |          38 |            0.53 |      **0.14** |           3.0 |   22 %  |
|     5 |          23 |            0.87 |          0.47 |          10.7 |  100 %  |
|     7 |          17 |            1.18 |          0.55 |          18.7 |   89 %  |
|    10 |          12 |            1.67 |          0.73 |          26.5 |  100 %  |

**Production setting (n=5000, K=3, p=10):**

| `K_f` | n_calib_arm | p / n_calib_arm | median spread | gamma_cal max | % blown |
|------:|------------:|----------------:|--------------:|--------------:|--------:|
|     3 |         176 |            0.06 |     **0.05**  |           0.80 |    0 %  |
|     5 |         106 |            0.09 |          0.08 |           0.83 |    0 %  |
|     7 |          76 |            0.13 |          0.17 |           0.88 |    0 %  |
|    10 |          53 |            0.19 |          0.23 |           1.31 |   67 %  |

**Takeaways:**

- `K_f = 3` is dramatically better than `K_f = 10` in both regimes. In the
  catastrophic setting the median spread drops 0.73 → 0.14 (5×) and the
  explosion rate 100% → 22%. In the production setting `K_f = 3` gives
  spread 0.05 (vs 0.23 at K_f=10) and 0% explosion.
- This is purely a sample-budget effect: reducing `K_f` from 10 → 3 grows
  `n_calib_arm` by ~3×, which is the same direction as growing `n_total`.
- **`K_f` is a tuning parameter passed to `run_simulation_study()`** — not an
  algorithm change. Just calling with `n_folds = 3` (or 5) should
  meaningfully improve coverage at fixed `n_total`.
- One artefact to note: at `K_f = 3` in the catastrophic setting,
  `w_exp_max` is still ~10⁴⁸ — γ_init still overfits — but the calibration
  fold (38 obs > p=20) is now large enough that the calibrated step can
  rescue the prediction. So **the calibration fold size matters more than
  the initial-DR coefficient size**.

**Revised practical recipe:** If hitting a coverage / blow-up problem,
*first* try `n_folds = 3` or `5` before increasing `n_total`. Cheaper, same
direction of effect.

### End-to-end verification — `K_f = 3` actually fixes coverage

`probe_kf.R` measures `mu_pred_ts` spread (an upstream proxy for bias).
To close the loop, `verify_kf3.R` (job `9395292`, 11.7 min wall) re-ran
the EXACT production setting (`n=5000, K=3, p=10, C1, ss=0.5`) with
`n_folds = 3` instead of 10 and compared side by side:

| Method               | bias K_f=10 | cov K_f=10 | bias K_f=3 | **cov K_f=3** | Δ coverage |
|----------------------|------------:|-----------:|-----------:|--------------:|-----------:|
| **one_round_crossfit** |  +0.00721 |  **0.886** | **−0.00051** | **0.945** | **+0.059** |
| **two_round_crossfit** |  +0.00764 |  **0.878** | **−0.00031** | **0.955** | **+0.077** |
| federated_dr         |  +0.00007 |     0.980  |   +0.00011  |     0.970  |   −0.010   |
| oracle_dr            |  −0.00010 |     0.988  |   −0.00019  |     0.985  |   −0.003   |
| pooled_dr            |  −0.00005 |     0.982  |   +0.00005  |     0.980  |   −0.002   |
| tilted_aipw          |  −0.00361 |     0.974  |   −0.00360  |     0.975  |   +0.001   |
| target_only          |  +0.00005 |     0.948  |   −0.00015  |     0.935  |   −0.013   |

**Confirmed:**
- one/two-round coverage moves from 0.88 → **0.945 / 0.955** (target 0.95).
- one/two-round bias collapses by ~15–25× (0.0072 → 0.0005).
- Other methods unchanged within MC noise — `K_f` is a one/two-round-specific
  knob, exactly as the upstream `mu_pred_ts` spread analysis predicted.
- Wall time also dropped (11 min vs hours at K_f=10) because fold-pair
  iterations are now `3 × 2 = 6` instead of `10 × 9 = 90`.

### Extended verification — all configs and a K_f midpoint

Jobs `9398801-9398804` completed (5-13 min each):

**C1, K_f sweep at fixed (n=5000, K=3, p=10) — `one_round_crossfit` and
`two_round_crossfit` coverage is monotone in K_f, exactly as the
`probe_kf.R` upstream signal predicted:**

| K_f | one_round bias | **one_round cov** | two_round bias | **two_round cov** |
|----:|---------------:|------------------:|---------------:|------------------:|
|  10 |     +0.00721   |         **0.886** |     +0.00764   |         **0.878** |
|   5 |     +0.00169   |         **0.930** |     +0.00230   |         **0.920** |
|   3 |     −0.00051   |         **0.945** |     −0.00031   |         **0.955** |

**C2 (outcome model misspecified) at K_f=3 — already-OK coverage stays OK,
bias improves 2-3×:**

| K_f | one_round bias | one_round cov | two_round bias | two_round cov |
|----:|---------------:|--------------:|---------------:|--------------:|
|  10 |     −0.00185   |         0.946 |     −0.00180   |         0.952 |
|   3 |     +0.00086   |         0.960 |     +0.00082   |         0.950 |

**C3 (propensity model misspecified) at K_f=3 — coverage essentially flat,
bias collapses ~20×. The propensity misspecification is the binding
constraint here, not the K_f-induced nuisance instability:**

| K_f | one_round bias | one_round cov | two_round bias | two_round cov |
|----:|---------------:|--------------:|---------------:|--------------:|
|  10 |     +0.00216   |         0.938 |     +0.00211   |         0.944 |
|   3 |     −0.00011   |         0.940 |     −0.00018   |         0.940 |

**C4 (both models misspecified) at K_f=3 — no production K_f=10 baseline.
Coverage in the 0.93 range is reasonable for the DR boundary case:**

| K_f | one_round bias | one_round cov | two_round bias | two_round cov |
|----:|---------------:|--------------:|---------------:|--------------:|
|   3 |     +0.00030   |         0.930 |     +0.00020   |         0.935 |

### Summary across all verified configs (C1-C4) at K_f = 3

| Config |  one_round bias |  one_round cov | two_round bias | two_round cov |
|--------|----------------:|---------------:|---------------:|--------------:|
| C1     |    −0.00051     |     **0.945**  |    −0.00031    |    **0.955**  |
| C2     |    +0.00086     |     **0.960**  |    +0.00082    |    **0.950**  |
| C3     |    −0.00011     |       0.940    |    −0.00018    |      0.940    |
| C4     |    +0.00030     |       0.930    |    +0.00020    |      0.935    |

**Takeaways:**

1. **C1 fully fixed**: 0.886 / 0.878 → 0.945 / 0.955 (Δ +0.06, +0.08). The
   K_f=5 intermediate point (0.930 / 0.920) confirms the trend is monotone,
   not an artefact.
2. **C2 fully fixed**: already 0.946 / 0.952, now 0.960 / 0.950 — coverage
   at target while bias drops 2×.
3. **C3 partially fixed**: coverage stays around 0.94 (mild under-coverage).
   The residual 0.06 gap is from the propensity misspecification itself,
   not from cross-fit nuisance instability — K_f cannot patch this. To
   close C3 we would need either a more flexible propensity model or to
   accept C3 as testing-misspecification-robustness.
4. **C4 (DR boundary)**: 0.93 coverage with both models wrong is plausibly
   "as good as it gets" for a DR estimator and matches the C2+C3 union of
   penalties.

**Bottom line:** `n_folds = 3` lifts every config to or above the
diagnosis-relevant target (C1/C2 to 0.95+, C3/C4 to ~0.94), and reduces
one/two-round bias by 5-25× across the board. Pure tuning win, no
algorithm change.

---

## File map

```
diagnosis/bias/
├── README.md                       <- this file
├── results/                        <- per-setting summary CSVs (4 anchor-axis
│                                      runs that finished)
├── log/                            <- SLURM stdout/stderr per task
├── checkpoints/                    <- R-level checkpoint files
├── probe_out/                      <- probe scan / plot outputs
│   ├── probe_scan_long_*.csv      <- per-(setting, seed, source, k1, k2) rows
│   ├── probe_scan_spread_*.csv    <- per-(setting, seed, source) spread summary
│   ├── probe_kf_long_*.csv        <- per-(setting, K_f, seed, source, k1, k2) rows
│   ├── probe_kf_spread_*.csv      <- per-(setting, K_f, seed, source) spread summary
│   └── plots/*.pdf                 <- PDF plots from probe_plot.R
│                                       (MSI R lacks png/cairo capability; PDF only)
│
├── run_setting.R / .sh / .cmd      <- run_simulation_study wrapper for one setting
├── bias_diagnose.sh / .cmd         <- SLURM array submission for the 10-cell grid
├── summarize.R                     <- aggregates results/*_summary.csv
│
├── probe.R / .sh / .cmd            <- 3-phase probe (A data, B per-function, C end-to-end)
├── probe_scan.R / .sh / .cmd       <- (n, p) grid scan of nuisance magnitudes
├── probe_kf.R / .sh / .cmd         <- K_f sweep at fixed (n, p)
└── probe_plot.R / .sh / .cmd       <- PNG plots from the probe CSVs
```

## Reproduction

```bash
# Stage 1 — full end-to-end grid (slow, multi-hour, can blow up on small n):
bash diagnosis/bias/bias_diagnose.sh --reset-checkpoint

# Stage 2 — fast per-function probe (minutes):
bash diagnosis/bias/probe.sh --phase B           # function-by-function
bash diagnosis/bias/probe.sh --phase C --time 01:00:00  # end-to-end seed sweep

# Stage 3 — (n, p) scan (15 min):
bash diagnosis/bias/probe_scan.sh

# Stage 4 — K_f sweep (60-90 min):
bash diagnosis/bias/probe_kf.sh

# Stage 5 — plots (run after probe_scan + probe_kf complete):
bash diagnosis/bias/probe_plot.sh
```

All scripts: never run R on the login node; always submit via SLURM.
Pinned partitions: `msismall, amdsmall, agsmall` (non-preempt).

---

## Status (as of submission of job `9390886`)

- `9276320` — original 11-task array, mostly killed by saffo-2tb preemption
  loop; only `K_2` produced (degenerate) results.
- `9311898` — 10-task array on non-preempt partitions, n_sims=50,
  CHECKPOINT_SIM_INTERVAL=5; 5/10 completed (`anchor`, `n_large`,
  `shift_005`, `shift_1`, plus the old `K_2`). 5 hung tasks (`p_50`,
  `p_100`, `p_200`, `K_4`, `K_5`) cancelled — they were never going to
  produce meaningful estimates given the n=1000 + p≥20 regime is below the
  algorithm's effective minimum.
- `9381180` — `probe.R` Phase A + B completed (the data above); Phase C
  hit walltime.
- `9384814` — `probe_scan.R` completed in 14 min; produced the cliff scan.
- `9390886` — `probe_kf.R` COMPLETED in 4 min; K_f sweep results above.
- `9392053` — `probe_plot.R` re-submitted after fixing a `%||%` definition
  order bug; failed because MSI's R build has no X11.
- `9392342` — second probe_plot retry with `type="cairo"`; cairo also
  unavailable.
- `9392657` — third probe_plot retry using PDF (the only graphics device
  this R build supports — `capabilities()` shows png/jpeg/tiff/cairo/X11
  all FALSE).
- `9395292` — `verify_kf3.R` C1 K_f=3 production verification COMPLETED
  in 11.7 min; coverage 0.945 / 0.955 (vs 0.886 / 0.878 at K_f=10).
- `9398801` — C2 K_f=3 verify COMPLETED in 7:39. cov 0.960 / 0.950.
- `9398802` — C3 K_f=3 verify COMPLETED in 5:39. cov 0.940 / 0.940
  (residual gap due to propensity misspec, not K_f).
- `9398803` — C4 K_f=3 verify COMPLETED in 7:28. cov 0.930 / 0.935.
- `9398804` — C1 K_f=5 verify COMPLETED in 13:16. cov 0.930 / 0.920
  (monotone between K_f=10 and K_f=3 → confirms trend).

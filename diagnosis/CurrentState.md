# Current state of <PROJECT NAME>

> **Single source of truth for what is true RIGHT NOW.**
> This file is rewritten (not appended) on every successful iteration.
> History lives in `HISTORY.md`.
>
> Last updated: <YYYY-MM-DD>
> Last updating HISTORY entry: <link>

---

## 1. Method default (1-line summary per knob)

- <core method or version, e.g., "TR-SP-HIP-v2 dirty Δ + χ² classifier">
- <hyperparameter defaults, e.g., "λ_Δ = 1.0, ω_y = auto = P">
- <any active config flags>

> Authoritative spec: `docs/method.tex` (Rule 13).  This section is a
> 1-line-per-knob restatement; the equations live in method.tex.

## 2. Headline numbers (multi-seed verified)

| Endpoint | Metric | Value | n_seeds | HISTORY ref |
|---|---|---|---|---|
| <e.g., sim 16-scenario> | <e.g., R²> | <e.g., 6W/6T/4L vs baselines> | <e.g., 5> | <link> |
| <e.g., TCGA binary> | <e.g., AUC> | <e.g., 0.792 ± 0.027> | <e.g., 5> | <link> |

## 3. Active issues (open iterations)

For each open iteration, list:

- **Symptom**: <one line>
- **Status**: PROPOSED / IN-FLIGHT / MEASURED awaiting decision
- **Stage**: 1 (diagnosis) / 2 (migration)  — Rule 16
- **HISTORY entry**: <link>
- **Blocker**: <e.g., waiting for sweep, waiting for review, waiting for input>

## 4. Parked for later (not v1 blockers)

- <topic 1>: <one line> (HISTORY: <link>)
- <topic 2>: <one line> (HISTORY: <link>)

## 5. Standing rules (or link to skill RULES.md)

Project-specific extensions (rule 20+):

- <none yet, OR list of project-specific promoted rules>

## 6. How to resume work (resume protocol)

To pick up this project after any interruption:

Total onboarding budget: ~10 minutes.

## 7. Project-specific resources

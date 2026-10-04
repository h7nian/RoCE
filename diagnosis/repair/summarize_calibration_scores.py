#!/usr/bin/env python3
"""Summarize candidate and eta-exchange diagnostics from verified saved scores."""
import argparse
from collections import defaultdict
import csv
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from review_current_calibration import paired_summary, require_close
from review_layer_comparison import metrics
from summarize_repeat_pilot import write_csv


COMPARISONS = (
    ("total", "full", "aggregate", "ablation", "aggregate"),
    ("candidate_at_full_eta", "full", "aggregate", "ablation_full_eta", "aggregate"),
    ("eta_on_ablation", "ablation_full_eta", "aggregate", "ablation", "aggregate"),
    ("eta_on_full", "full", "aggregate", "full_ablation_eta", "aggregate"),
    ("candidate_at_ablation_eta", "full_ablation_eta", "aggregate", "ablation", "aggregate"),
    ("source_s1", "full", "s1", "ablation", "s1"),
    ("source_s2", "full", "s2", "ablation", "s2"))


def summarize(root):
    root = pilot.scratch_path(root)
    if not (root / "COMPLETE").is_file():
        raise ValueError("Require completed score reconstruction")
    with (root / "score_estimates.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    groups = defaultdict(dict)
    expected = {(variant, "aggregate", quantity)
                for variant in ("full", "ablation", "ablation_full_eta", "full_ablation_eta")
                for quantity in ("mu1", "mu0", "tate")}
    expected |= {(variant, candidate, quantity)
                 for variant in ("full", "ablation") for candidate in ("target_anchor", "s1", "s2")
                 for quantity in ("mu1", "mu0", "tate")}
    dimensions = {int(row["p"]) for row in rows}
    if len(dimensions) != 1 or not dimensions.issubset({100, 200}):
        raise ValueError("One high-dimensional panel is required")
    for row in rows:
        if int(row["K"]) != 2 or float(row["rho"]) != 0:
            raise ValueError("Unexpected calibration diagnostic setting")
        key = row["config"], row["variant"], row["candidate"], row["quantity"]
        seed = int(row["sim_id"])
        if seed in groups[key]:
            raise ValueError("Duplicate score row")
        groups[key][seed] = row
    if set(groups) != {(config,) + key for config in ("C1", "C2", "C3") for key in expected}:
        raise ValueError("Missing or extra score diagnostic group")
    for indexed in groups.values():
        if set(indexed) != set(range(1, 201)):
            raise ValueError("Each score diagnostic requires all200 prescribed seeds")
    for config in ("C1", "C2", "C3"):
        for seed in range(1, 201):
            hashes = {indexed[seed]["data_sha256"] for key, indexed in groups.items() if key[0] == config}
            if len(hashes) != 1:
                raise ValueError("Score diagnostic data hashes differ")
        for quantity in ("mu1", "mu0", "tate"):
            for seed in range(1, 201):
                for field in ("estimate", "truth", "variance", "variance_fixed_weights"):
                    require_close(groups[config, "full", "target_anchor", quantity][seed][field],
                                  groups[config, "ablation", "target_anchor", quantity][seed][field],
                                  "unchanged target anchor " + field)
    summaries, comparisons = [], []
    for (config, variant, candidate, quantity), indexed in sorted(groups.items()):
        values = [indexed[seed] for seed in range(1, 201)]
        summaries.append(dict(config=config, variant=variant, candidate=candidate, quantity=quantity,
                              **metrics(values, 200)))
    for config in ("C1", "C2", "C3"):
        for quantity in ("mu1", "mu0", "tate"):
            contrasts = {}
            for label, left_variant, left_candidate, right_variant, right_candidate in COMPARISONS:
                value = paired_summary(
                    [groups[config, left_variant, left_candidate, quantity][seed] for seed in range(1, 201)],
                    [groups[config, right_variant, right_candidate, quantity][seed] for seed in range(1, 201)])
                contrasts[label] = value
                comparisons.append(dict(config=config, quantity=quantity, comparison=label, **value))
            for first, second in (("candidate_at_full_eta", "eta_on_ablation"),
                                  ("eta_on_full", "candidate_at_ablation_eta")):
                require_close(contrasts[first]["mse_difference"] + contrasts[second]["mse_difference"],
                              contrasts["total"]["mse_difference"], "MSE path decomposition")
    output = root / "summary"
    output.mkdir(exist_ok=False)
    write_csv(output / "metrics.csv", summaries)
    write_csv(output / "paired_comparisons.csv", comparisons)
    shutil.copy2(__file__, output / Path(__file__).name)
    pilot.write_json(output / "checks.json", dict(paired_datasets=600, groups=len(groups),
        input_sha256=pilot.digest(root / "score_estimates.csv"),
        code_sha256=pilot.digest(Path(__file__)),
        target_anchor_invariants=True, mse_path_decompositions=True))
    lines = ["# Candidate and eta diagnostic", "",
        "All600 paired datasets have identical hashes and actual outer-fold assignments. Reconstructed own-weight estimates and both variance formulas match the committed estimator within1e-12; see ../checks.txt.",
        "Every cell has200 prescribed seeds. Negative paired MSE difference favors the first program. Weight-exchange rows are diagnostics using both fitted programs, not proposed new estimators.", "",
        "| Scenario | Contrast | Paired TATE MSE difference | MCSE | Difference / MCSE |", "|---|---|---:|---:|---:|"]
    for row in comparisons:
        if row["quantity"] == "tate":
            lines.append("|{config}|{comparison}|{mse_difference:.3g}|{mse_difference_mcse:.3g}|{mse_difference_z:.2f}|".format(**row))
    lines += ["", "The total MSE difference equals candidate_at_full_eta + eta_on_ablation, and also eta_on_full + candidate_at_ablation_eta. These are algebraic, path-dependent decompositions, not independent causal contributions; their MCSEs must not be added.",
        "Source_s1/source_s2 compare each individual assisted estimator before aggregation. The target anchor is unchanged in both programs for every seed, arm, estimate and score variance.",
        "The full score variance includes the existing eta-sensitivity calculation holding nuisance fits and active sets fixed. Fixed-weight candidate variances are separately labeled; neither is a full nuisance-refit variance.", ""]
    (output / "report.md").write_text("\n".join(lines))
    print("Verified600 paired score datasets; {} summaries and {} contrasts.".format(len(summaries),len(comparisons)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    summarize(parser.parse_args().root)

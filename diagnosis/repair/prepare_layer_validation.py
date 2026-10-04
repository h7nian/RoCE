#!/usr/bin/env python3
"""Prepare paired two/three-layer experiments with identical outer estimators."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import EXCLUDED_NODES, experiment_identity, first_seed_gate, read_campaign


def profiles(seed_start=401, n_repeats=50):
    if seed_start < 1 or n_repeats < 1:
        raise ValueError("Positive seed start and repeat count are required")
    return [dict(name="layers{}_p{}".format(layers, dimension), crossfit_layers=layers,
                 dimension=dimension, sources=[2, 4, 6], configs=["C1", "C2", "C3"],
                 rhos=[0], seed_start=seed_start, n_repeats=n_repeats,
                 max_in_flight=32 if dimension == 10 else 48)
            for dimension in (10, 100, 200) for layers in (2, 3)]


def prepare(output, check_root, pilot_root, seed_start=401, n_repeats=50):
    output, check_root, pilot_root = map(pilot.scratch_path, (output, check_root, pilot_root))
    for marker in ("TESTS_PASSED", "CROSSFIT_LAYER_CHECKS_PASSED"):
        if not (check_root / marker).is_file():
            raise ValueError("Missing checked-library marker: " + marker)
    all_profiles = profiles(seed_start, n_repeats)
    output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    workflow = output / "coordinator_workflow"
    workflow.mkdir()
    for name in (Path(__file__).name, "submit_repeat_pilot.py", "advance_validation_plan.py", "prepare_highdim_validation.py"):
        shutil.copy2(source / name, workflow / name)
    pilot.write_json(workflow / "manifest.json", {str(path): pilot.digest(path) for path in workflow.iterdir()})
    pilot.write_json(output / "profiles.json", all_profiles)
    initial_gate = dict(root=str(pilot_root / "layers2"), paired_root=str(pilot_root / "layers3"), tasks=[1])
    campaigns, registry, seen = [], [], set()
    for profile in all_profiles:
        root = output / profile["name"]
        layers, dimension = profile["crossfit_layers"], profile["dimension"]
        command = [sys.executable, "-B", str(source / "submit_repeat_pilot.py"), "prepare", str(root),
                   "--check-root", str(check_root), "--crossfit-layers", str(layers),
                   "--dimensions", str(dimension), "--configs"] + profile["configs"]
        command += ["--sources"] + [str(value) for value in profile["sources"]]
        command += ["--rhos", "0", "--seed-start", str(seed_start), "--n-repeats", str(n_repeats),
                    "--recipes", "score_derivative", "--aggregation-cutoff", "2", "--cpus", "2",
                    "--memory", "4G" if dimension == 10 else "12G",
                    "--time-limit", "00:30:00" if dimension == 10 else "02:00:00",
                    "--exclude-nodes", EXCLUDED_NODES]
        if dimension == 10:
            command += ["--completed-from", str(pilot_root / ("layers" + str(layers)))]
        subprocess.run(command, check=True, env=pilot.submission_environment())
        configuration, rows = read_campaign(root)
        for row in rows:
            identity = experiment_identity(row, configuration)
            if identity in seen:
                raise ValueError("Duplicated layer experiment: " + str(identity))
            seen.add(identity)
        gates = [initial_gate]
        if dimension != 10:
            first_root = output / "layers2_p10"
            _, first_rows = read_campaign(first_root)
            paired_gate = first_seed_gate(first_root, first_rows)
            paired_gate["paired_root"] = str(output / "layers3_p10")
            gates.append(paired_gate)
        campaigns.append(dict(name=root.name, root=str(root), max_in_flight=profile["max_in_flight"],
                              theory_role="fixed_K_crossfit_layer_ablation", gates=gates))
        registry.append(dict(root=str(root), planned=len(rows), crossfit_layers=layers,
                             dimension=dimension, origin="new"))
    # Pair on the entire manifest, not only the first seed used by release gates.
    for dimension in (10, 100, 200):
        _, two = read_campaign(output / "layers2_p{}".format(dimension))
        _, three = read_campaign(output / "layers3_p{}".format(dimension))
        if two != three:
            raise ValueError("Two/three-layer manifests are not paired")
    plan = dict(purpose="Compare eta, both potential-outcome means, TATE and variance with fixed K; no performance-based gates.",
                profile_set="layers", new_repeats=dict(method=sum(row["planned"] for row in registry)),
                max_new_workers=sum(row["max_in_flight"] for row in campaigns), campaigns=campaigns)
    pilot.write_json(output / "expansion_plan.json", plan)
    pilot.write_json(output / "campaign_registry.json", registry)
    pilot.write_json(output / "design_checks.json", dict(duplicate_experiments=0, paired_manifests=True,
        n_per_site=1000, n_folds=10, nuisance_grid=100, cutoff=2, rhos=[0],
        target_benchmarks="Same-repeat ordinary and calibrated target-only estimates from each worker",
        scope="Initial {}-pair variance comparison; report Monte Carlo uncertainty. Deviations and growing K remain separate experiments.".format(n_repeats)))
    print("Prepared {} paired layer tasks with a maximum of {} new workers".format(len(seen), plan["max_new_workers"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-root", required=True, type=Path)
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--seed-start", type=int, default=401)
    parser.add_argument("--n-repeats", type=int, default=50)
    args = parser.parse_args()
    prepare(args.output, args.check_root, args.pilot_root, args.seed_start, args.n_repeats)

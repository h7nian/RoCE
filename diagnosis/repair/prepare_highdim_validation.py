#!/usr/bin/env python3
"""Prepare a theory-motivated, nonduplicated high-dimensional validation plan.

This command prepares frozen campaigns and release gates. The existing
advance_validation_plan.py submits controllers after those gates pass.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot

EXCLUDED_NODES = "acl03,acl17,acl18,acl19,acl22,acl41,acl45,acl47,acl53,acl66,acl72,acl73,acl74,acl83,acl96,cn1114"
ABLATIONS = {
    "target_lasso": {"target_nuisance_method": "lasso"},
    "initial_validation": {"source_validation_method": "initial"},
    "two_round": {"protocol": "two_round"},
}


def profiles(profile_set="highdim"):
    result = []

    def add(name, dimensions, rhos, role, maximum, baseline_maximum=0, **options):
        profile = dict(name=name, dimensions=dimensions, rhos=rhos, role=role,
                       configs=["C1", "C2", "C3"], seed_start=1, n_repeats=50,
                       sources=[2], n_deviated_sites=1, memory="8G",
                       protocol="one_round", target_nuisance_method="hou_calibrated",
                       source_nuisance_method="calibrated", source_validation_method="calibrated", covariate_shift=1,
                       source_treatment_scale=1, max_in_flight=maximum,
                       baseline_max_in_flight=baseline_maximum)
        profile.update(options)
        result.append(profile)

    if profile_set == "calibration":
        for target in ("hou_calibrated", "lasso"):
            name = "standard_source" + ("_target_lasso" if target == "lasso" else "")
            previous = None
            for dimension, repeats, maximum in ((10, 1, 3), (100, 50, 32), (200, 50, 32)):
                current = name + "_p" + str(dimension)
                options = dict(source_nuisance_method="standard", source_validation_method="complete",
                               target_nuisance_method=target, n_repeats=repeats)
                if previous is not None:
                    options["parent"] = previous
                add(current, [dimension], [0], "final_score_calibration_ablation", maximum, **options)
                previous = current
        return result
    if profile_set == "weak_fraction":
        for sources, deviated in ((4, 2), (6, 3)):
            add("half_invalid_k{}_p100_rho_half".format(sources), [100], [.5],
                "weak_fixed_invalid_fraction", 48, 24, sources=[sources],
                n_deviated_sites=deviated, memory="12G" if sources == 6 else "8G",
                existing_parent="source_count_validation_v1/half_invalid_k{}_p100".format(sources))
        return result
    if profile_set == "source_count":
        add("source_count_p100", [100], [0, .5, 1, 2], "fixed_invalid_count", 96, 48,
            sources=[4, 6], memory="12G")
        add("source_count_p200", [200], [0, .5, 1, 2], "fixed_invalid_count", 96, 48,
            sources=[4, 6], memory="12G")
        add("half_invalid_k4_p100", [100], [1, 2], "fixed_invalid_fraction", 48, 24,
            sources=[4], n_deviated_sites=2)
        add("half_invalid_k6_p100", [100], [1, 2], "fixed_invalid_fraction", 48, 24,
            sources=[6], n_deviated_sites=3, memory="12G")
        return result
    if profile_set != "highdim":
        raise ValueError("Unknown profile set: " + profile_set)
    add("outcome_shift_p200", [200], [.5, 1, 2], "deviation_validation", 128, 64)
    add("strong_deviation_p100_seeds51_200", [100], [1, 2], "primary", 128, 64,
        seed_start=51, n_repeats=150, existing_parent="extended_validation_v1/outcome_shift")
    add("strong_deviation_p200_seeds51_200", [200], [1, 2], "primary", 128, 64,
        seed_start=51, n_repeats=150, parent="outcome_shift_p200")
    add("covariate_shift_half_p200", [200], [0], "transport_sensitivity", 48, 16,
        covariate_shift=.5)
    add("covariate_shift_double_p200", [200], [0], "transport_sensitivity", 48, 16,
        covariate_shift=2)
    for name, options in ABLATIONS.items():
        add(name + "_p200", [200], [0, 1], "mechanism_ablation", 40,
            parent="outcome_shift_p200", **options)
        add(name + "_rho2", [100, 200], [2], "mechanism_ablation", 40,
            parent="outcome_shift_p200", **options)
    add("both_models_misspecified", [100, 200], [0, 1, 2], "outside_nuisance_guarantee", 32, 16,
        configs=["C4"])
    add("identical_populations", [100, 200], [0], "transport_negative_control", 16, 16,
        configs=["C1"], covariate_shift=0, source_treatment_scale=0)
    return result


def rows_for_profile(profile):
    return pilot.task_rows(profile["configs"], profile["sources"], profile["rhos"],
        range(profile["seed_start"], profile["seed_start"] + profile["n_repeats"]),
        profile["protocol"], "treated_arm", profile["dimensions"])


def data_identity(row, configuration):
    control = configuration.get("dgp_control") or {}
    declared_count = int(configuration.get("n_deviated_sites", 1))
    rho = float(row["rho"]) if declared_count else 0.0
    active_count = declared_count if rho else 0
    return (str(row["config"]), int(row["K"]), int(row["p"]), rho, active_count,
            int(row["sim_id"]), str(row["deviation_mechanism"]),
            configuration.get("dgp_type", "bounded"), configuration.get("n_per_site", 1000),
            float(control.get("covariate_shift", 1)), float(control.get("source_treatment_scale", 1)))


def experiment_identity(row, configuration):
    data = data_identity(row, configuration)
    shared = (configuration.get("n_folds", 10), configuration.get("nlambda", 100),
              configuration.get("nuisance_tol", 1e-10))
    if configuration.get("baseline_methods"):
        return ("baseline", data, shared, tuple(sorted(configuration["baseline_methods"])),
                configuration.get("baseline_variance_method", "bootstrap"),
                configuration.get("n_bootstrap", 5000))
    validation = configuration.get("source_validation_method", "calibrated")
    if configuration.get("crossfit_layers") == 2:
        validation = "outer_fit"
    elif configuration.get("crossfit_layers") == 3 and validation == "outer_fit":
        raise ValueError("outer_fit weight learning requires crossfit_layers=2")
    return ("method", data, shared, str(row["protocol"]),
            configuration.get("target_nuisance_method", "hou_calibrated"),
            configuration.get("source_nuisance_method", "calibrated"),
            "complete" if validation in ("complete", "calibrated") else validation,
            tuple(configuration.get("recipes", ["score_derivative"])),
            configuration.get("source_radius", 12), configuration.get("target_score_radius", 12),
            configuration.get("aggregation_lambda", .5))


def read_campaign(root):
    configuration = json.loads((root / "configuration.json").read_text())
    with (root / "manifest.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    return configuration, rows


def first_seed_gate(root, rows):
    first = min(int(row["sim_id"]) for row in rows)
    return dict(root=str(root), tasks=[int(row["task_id"]) for row in rows if int(row["sim_id"]) == first])


def prepare(base, output, profile_set="highdim", method_check_root=None):
    base, output = pilot.scratch_path(base), pilot.scratch_path(output)
    if profile_set == "calibration" and method_check_root is None:
        raise ValueError("The calibration ablation requires its checked installed library via --check-root")
    method_check_root = pilot.scratch_path(method_check_root or base / "checks_v1")
    existing = [base / "validation_highdim_mc200_cutoff2_v1", base / "baseline_mc200_v1"]
    existing += sorted(root for root in (base / "extended_validation_v1").iterdir()
                       if (root / "manifest.csv").exists() and not root.name.endswith("_check"))
    if profile_set in ("source_count", "weak_fraction"):
        earlier = json.loads((base / "highdim_validation_v2/campaign_registry.json").read_text())
        existing += [Path(entry["root"]) for entry in earlier if entry["origin"] == "new"]
    if profile_set == "weak_fraction":
        earlier = json.loads((base / "source_count_validation_v1/campaign_registry.json").read_text())
        existing += [Path(entry["root"]) for entry in earlier if entry["origin"] == "new"]
    seen, baseline_data, registry = {}, set(), []
    for root in existing:
        configuration = pilot.verify(root)
        _, rows = read_campaign(root)
        for row in rows:
            identity = experiment_identity(row, configuration)
            if identity in seen:
                raise ValueError("Existing production campaigns overlap: " + str((root, seen[identity])))
            seen[identity] = str(root)
            if configuration.get("baseline_methods"):
                baseline_data.add(data_identity(row, configuration))
        registry.append(dict(root=str(root), kind="baseline" if configuration.get("baseline_methods") else "method",
                             origin="existing", planned=len(rows)))

    output.mkdir(parents=True, exist_ok=False)
    workflow = output / "coordinator_workflow"
    workflow.mkdir()
    source = Path(__file__).resolve().parent
    for name in (Path(__file__).name, "submit_repeat_pilot.py", "advance_validation_plan.py"):
        shutil.copy2(source / name, workflow / name)
    pilot.write_json(workflow / "manifest.json", {str(path): pilot.digest(path) for path in workflow.iterdir()})

    common_gates = [dict(root=str(base / "validation_highdim_mc200_cutoff2_v1"), tasks=list(range(1, 7))),
                    dict(root=str(base / "baseline_pilot_v2"), tasks=list(range(1, 10)))]
    if profile_set in ("source_count", "weak_fraction"):
        checks_root = base / "source_count_checks_v1"
        expectations = json.loads((checks_root / "probe_data_expectations.json").read_text())
        for name, expected in expectations.items():
            common_gates.append(dict(root=str(checks_root / name), tasks=[1],
                                     expected_data_hashes={"1": expected}))
    campaign_plan, method_data = [], set()
    all_profiles = profiles(profile_set)
    pilot.write_json(output / "profiles.json", all_profiles)
    for profile in all_profiles:
        name = profile["name"]
        for baseline in (False, True):
            maximum = profile["baseline_max_in_flight"] if baseline else profile["max_in_flight"]
            if maximum == 0:
                continue
            root = output / (name + ("_baselines" if baseline else ""))
            command = [sys.executable, "-B", str(source / "submit_repeat_pilot.py"), "prepare", str(root),
                       "--check-root", str(base / "baseline_checks_v1" if baseline else method_check_root),
                       "--dimensions"] + [str(value) for value in profile["dimensions"]]
            command += ["--configs"] + profile["configs"] + ["--sources"]
            command += [str(value) for value in profile["sources"]] + ["--rhos"]
            command += [str(value) for value in profile["rhos"]]
            command += ["--seed-start", str(profile["seed_start"]), "--n-repeats", str(profile["n_repeats"]),
                        "--aggregation-cutoff", "2", "--covariate-shift", str(profile["covariate_shift"]),
                        "--source-treatment-scale", str(profile["source_treatment_scale"]),
                        "--n-deviated-sites", str(profile["n_deviated_sites"]),
                        "--protocol", profile["protocol"], "--cpus", "2", "--exclude-nodes", EXCLUDED_NODES,
                        "--memory", "4G" if baseline else profile["memory"],
                        "--time-limit", "01:00:00" if baseline else "02:00:00"]
            if baseline:
                command += ["--baseline-methods"] + list(pilot.BASELINE_METHODS)
                command += ["--baseline-reference-roots", str(output / name)]
            else:
                command += ["--recipes", "score_derivative", "--target-nuisance-method", profile["target_nuisance_method"],
                            "--source-nuisance-method", profile["source_nuisance_method"],
                            "--source-validation-method", profile["source_validation_method"]]
            subprocess.run(command, check=True, env=pilot.submission_environment())
            configuration, rows = read_campaign(root)
            for row in rows:
                identity = experiment_identity(row, configuration)
                if identity in seen:
                    raise ValueError("New campaign duplicates an existing experiment: " + str((root, seen[identity])))
                seen[identity] = str(root)
                (baseline_data if baseline else method_data).add(data_identity(row, configuration))
            gates = list(common_gates)
            if profile.get("parent"):
                parent_root = output / profile["parent"]
                _, parent_rows = read_campaign(parent_root)
                gates.append(first_seed_gate(parent_root, parent_rows))
            if profile.get("existing_parent"):
                parent_root = base / profile["existing_parent"]
                _, parent_rows = read_campaign(parent_root)
                gates.append(first_seed_gate(parent_root, parent_rows))
            if baseline:
                gates.append(first_seed_gate(output / name, rows))
            campaign_plan.append(dict(name=root.name, root=str(root), max_in_flight=maximum,
                                      theory_role=profile["role"], gates=gates))
            registry.append(dict(root=str(root), kind="baseline" if baseline else "method",
                                 origin="new", planned=len(rows)))

    if not method_data.issubset(baseline_data):
        raise ValueError("Some method data settings lack a planned or existing matched baseline")
    totals = Counter()
    for entry in registry:
        if entry["origin"] == "new":
            totals[entry["kind"]] += entry["planned"]
    concurrency = sum(campaign["max_in_flight"] for campaign in campaign_plan)
    expected_concurrency = {"highdim": 1024, "source_count": 432, "weak_fraction": 144, "calibration": 134}[profile_set]
    if concurrency != expected_concurrency:
        raise ValueError("Unexpected aggregate worker-concurrency budget")
    plan = dict(purpose="User-authorized high-dimensional theory checks and boundary experiments; no performance-based release gates.",
                profile_set=profile_set, new_repeats=dict(totals), max_new_workers=concurrency, campaigns=campaign_plan)
    pilot.write_json(output / "expansion_plan.json", plan)
    pilot.write_json(output / "campaign_registry.json", registry)
    pilot.write_json(output / "design_checks.json", dict(
        duplicate_experiments=0, method_settings_with_matched_baselines=len(method_data),
        all_new_method_data_have_baselines=True, new_repeats=dict(totals), max_new_workers=concurrency))
    print("Prepared {} campaigns: {}; maximum {} simultaneous new workers".format(
        len(campaign_plan), dict(totals), concurrency))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--profile-set", choices=["highdim", "source_count", "weak_fraction", "calibration"], default="highdim")
    parser.add_argument("--check-root", type=Path, help="Checked method library; required for the calibration ablation")
    arguments = parser.parse_args()
    prepare(arguments.run_root, arguments.output, arguments.profile_set, arguments.check_root)

#!/usr/bin/env python3
"""Review paired layer experiments, preserving failures and Monte Carlo uncertainty."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from summarize_repeat_pilot import write_csv, wilson_interval
from advance_validation_plan import validate_layer_pair_configurations

CELL = ("p", "K", "config", "rho", "protocol", "deviation_mechanism")
MODES = ("common_tate", "separate_arms", "joint_tate")
QUANTITIES = ("mu1", "mu0", "tate")


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def metrics(rows, planned):
    result = dict(n_planned=planned, n_success=len(rows), n_unavailable=planned-len(rows))
    if not rows:
        return result
    truth = [float(row["truth"]) for row in rows]
    if max(truth) - min(truth) > 1e-12:
        raise ValueError("Population truth differs within a layer cell")
    errors = [float(row["estimate"]) - value for row, value in zip(rows, truth)]
    empirical = statistics.variance(errors) if len(errors) > 1 else None
    result.update(bias=statistics.mean(errors), rmse=math.sqrt(statistics.mean(value*value for value in errors)),
                  empirical_variance=empirical,
                  bias_mcse=math.sqrt(empirical/len(rows)) if empirical is not None else None)
    for field, suffix in (("variance", ""), ("variance_fixed_weights", "_fixed_weights")):
        variances = [float(row[field]) for row in rows]
        if any(not math.isfinite(value) or value <= 0 for value in variances):
            raise ValueError("Invalid reported variance")
        covered = sum(abs(error) <= 1.959963984540054*math.sqrt(variance)
                      for error, variance in zip(errors, variances))
        lower, upper = wilson_interval(covered, len(rows))
        result.update({"mean_variance"+suffix: statistics.mean(variances),
                       "variance_ratio"+suffix: statistics.mean(variances)/empirical if empirical else None,
                       "coverage"+suffix: covered/len(rows),
                       "coverage_lower"+suffix: lower, "coverage_upper"+suffix: upper,
                       "mean_ci_width"+suffix: 2*1.959963984540054*statistics.mean(math.sqrt(value) for value in variances)})
    return result


def write_union_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    write_csv(path, [{key: row.get(key) for key in fields} for row in rows])


def paired_variance_comparison(first, second, reported_first, reported_second):
    count = len(first)
    if count == 0 or any(len(values) != count for values in (second, reported_first, reported_second)):
        raise ValueError("Variance comparison requires nonempty paired samples")
    reported = [x-y for x,y in zip(reported_first,reported_second)]
    mean_reported = statistics.mean(reported)
    empirical = statistics.variance(first)-statistics.variance(second) if count > 1 else None
    jackknife_se = None
    if count > 2:
        first_mean, second_mean = statistics.mean(first), statistics.mean(second)
        first_sse = sum((value-first_mean)**2 for value in first)
        second_sse = sum((value-second_mean)**2 for value in second)
        reported_sum = sum(reported)
        leave_one_out = []
        for x,y,variance in zip(first,second,reported):
            empirical_without = (first_sse-second_sse-count/(count-1)*((x-first_mean)**2-(y-second_mean)**2))/(count-2)
            reported_without = (reported_sum-variance)/(count-1)
            leave_one_out.append(reported_without-empirical_without)
        center = statistics.mean(leave_one_out)
        jackknife_se = math.sqrt((count-1)/count*sum((value-center)**2 for value in leave_one_out))
    return dict(empirical_variance_difference_two_minus_three=empirical,
                mean_reported_variance_difference_two_minus_three=mean_reported,
                variance_difference_discrepancy=mean_reported-empirical if empirical is not None else None,
                variance_discrepancy_jackknife_se=jackknife_se)


def covariance_summary(records):
    output = []
    for (cell, recipe, mode, layers), values in sorted(records.items()):
        first = [float(value["mu1"]["estimate"]) for value in values]
        second = [float(value["mu0"]["estimate"]) for value in values]
        first_mean, second_mean = statistics.mean(first), statistics.mean(second)
        empirical = sum((x-first_mean)*(y-second_mean) for x,y in zip(first,second))/(len(first)-1) if len(first)>1 else None
        output.append(dict(zip(CELL, cell), recipe=recipe, aggregation_mode=mode, crossfit_layers=layers,
            n_success=len(values), empirical_covariance=empirical,
            mean_covariance=statistics.mean(float(value["mu1"]["covariance_mu1_mu0"]) for value in values),
            mean_covariance_fixed_weights=statistics.mean(float(value["mu1"]["covariance_mu1_mu0_fixed_weights"]) for value in values)))
    return output


def covariance_accuracy(first, second, reported):
    """Compare mean reported covariance with Monte Carlo covariance.

    Delete-one jackknife uncertainty keeps the two estimates and their reported
    covariance paired within each independent simulation repeat.
    """
    count = len(first)
    if count < 2 or len(second) != count or len(reported) != count:
        raise ValueError("Covariance accuracy requires at least two paired repeats")
    first_mean, second_mean = statistics.mean(first), statistics.mean(second)
    centered = [(x-first_mean)*(y-second_mean) for x,y in zip(first,second)]
    cross_sum, reported_sum = sum(centered), sum(reported)
    empirical, estimated = cross_sum/(count-1), reported_sum/count
    standard_error = None
    if count > 2:
        deleted = [(reported_sum-value)/(count-1) -
                   (cross_sum-count*product/(count-1))/(count-2)
                   for product,value in zip(centered,reported)]
        center = statistics.mean(deleted)
        standard_error = math.sqrt((count-1)/count*sum((value-center)**2 for value in deleted))
    return dict(empirical_component=empirical, mean_reported_component=estimated,
                reported_minus_empirical=estimated-empirical,
                discrepancy_jackknife_se=standard_error)


def variance_component_summary(records):
    output = []
    for (cell, recipe, mode, layers), values in sorted(records.items()):
        if len(values) < 2:
            continue
        first = [float(value["mu1"]["estimate"]) for value in values]
        second = [float(value["mu0"]["estimate"]) for value in values]
        contrast = [x-y for x,y in zip(first,second)]
        for variance, covariance, inference in (
                ("variance", "covariance_mu1_mu0", "eta_sensitivity"),
                ("variance_fixed_weights", "covariance_mu1_mu0_fixed_weights", "fixed_weights")):
            first_reported = [float(value["mu1"][variance]) for value in values]
            second_reported = [float(value["mu0"][variance]) for value in values]
            cross_reported = [float(value["mu1"][covariance]) for value in values]
            components = (("mu1", first, first, first_reported),
                          ("mu0", second, second, second_reported),
                          ("minus_twice_covariance", first, [-2*x for x in second], [-2*x for x in cross_reported]),
                          ("tate", contrast, contrast, [x+y-2*z for x,y,z in zip(first_reported,second_reported,cross_reported)]))
            for component, x, y, reported in components:
                output.append(dict(zip(CELL,cell), recipe=recipe, aggregation_mode=mode,
                    crossfit_layers=layers, n_pairs=len(values), inference=inference,
                    component=component, **covariance_accuracy(x,y,reported)))
    return output


def review(plan_root, output):
    plan_root, output = map(pilot.scratch_path, (plan_root, output))
    output.mkdir(parents=True, exist_ok=False)
    plan = json.loads((plan_root / "expansion_plan.json").read_text())
    groups, planned, pair_records, covariance_records, weight_records = defaultdict(list), defaultdict(int), {}, defaultdict(list), defaultdict(list)
    statuses, runtime = [], defaultdict(list)
    for campaign in plan["campaigns"]:
        root = pilot.scratch_path(campaign["root"])
        config = pilot.verify(root)
        layers = config["crossfit_layers"]
        for task in read_csv(root / "manifest.csv"):
            cell = tuple(task[key] for key in CELL)
            directory = root / "tasks" / task["task_id"]
            complete = (directory / "COMPLETE").is_file()
            state = (directory / "status.txt").read_text().strip() if (directory / "status.txt").is_file() else "NOT_STARTED"
            statuses.append(dict(task, crossfit_layers=layers, campaign=root.name, state=state, committed=complete))
            for recipe in config["recipes"]:
                for mode in MODES:
                    for quantity in QUANTITIES:
                        planned[(cell, recipe, mode, quantity, layers)] += 1
            if not complete:
                continue
            if state != "COMPLETE":
                raise ValueError("Inconsistent completion state: " + str(directory))
            rows = read_csv(directory / "layer_estimates.csv")
            expected = {(recipe, mode, quantity) for recipe in config["recipes"] for mode in MODES for quantity in QUANTITIES}
            if len(rows) != len(expected) or {(row["recipe"], row["aggregation_mode"], row["quantity"]) for row in rows} != expected:
                raise ValueError("Incomplete or duplicated layer outputs: " + str(directory))
            data_hash = (directory / "data_sha256.txt").read_text().strip()
            for row in rows:
                if int(row["crossfit_layers"]) != layers or any(not math.isfinite(float(row[key])) for key in
                        ("estimate", "truth", "variance", "variance_fixed_weights")):
                    raise ValueError("Invalid layer output: " + str(directory))
                key = (cell, row["recipe"], row["aggregation_mode"], row["quantity"], layers)
                groups[key].append(row)
                pair_key = key[:-1] + (task["sim_id"],)
                if layers in pair_records.setdefault(pair_key, {}):
                    raise ValueError("Duplicate layer/seed record: " + str(pair_key))
                pair_records[pair_key][layers] = dict(row, data_hash=data_hash,
                    configuration=config,
                    outer_hash=(directory / (row["recipe"]+"_outer_fit_sha256.txt")).read_text().strip())
            for recipe in config["recipes"]:
                for mode in MODES:
                    values = {row["quantity"]: row for row in rows if row["recipe"] == recipe and row["aggregation_mode"] == mode}
                    for field, covariance in (("variance", "covariance_mu1_mu0"),
                                              ("variance_fixed_weights", "covariance_mu1_mu0_fixed_weights")):
                        reconstructed = float(values["mu1"][field])+float(values["mu0"][field])-2*float(values["mu1"][covariance])
                        if abs(reconstructed-float(values["tate"][field])) > 1e-12:
                            raise ValueError("Arm covariance does not reproduce TATE variance")
                    if abs(float(values["mu1"]["estimate"])-float(values["mu0"]["estimate"])-float(values["tate"]["estimate"])) > 1e-12:
                        raise ValueError("Arm means do not reproduce TATE")
                    covariance_records[(cell, recipe, mode, layers)].append(values)
            for row in read_csv(directory / "layer_weights.csv"):
                key = (cell, row["recipe"], row["aggregation_mode"], row["coordinate"], layers, task["sim_id"])
                weight = float(row["weight"])
                if not math.isfinite(weight):
                    raise ValueError("Nonfinite aggregation weight")
                weight_records[key].append(weight)
            result_rows = read_csv(directory / "results.csv")
            runtime[(cell, layers)].append(float(result_rows[0]["simulation_elapsed_seconds"]))
    summaries = []
    for key, count in sorted(planned.items()):
        cell, recipe, mode, quantity, layers = key
        summaries.append(dict(zip(CELL, cell), recipe=recipe, aggregation_mode=mode, quantity=quantity,
                              crossfit_layers=layers, **metrics(groups[key], count)))
    paired = defaultdict(list)
    paired_groups = defaultdict(list)
    paired_covariance_records = defaultdict(list)
    for key, pair in pair_records.items():
        if set(pair) != {2, 3}:
            continue
        first, second = pair[2], pair[3]
        validate_layer_pair_configurations(first["configuration"], second["configuration"])
        if first["data_hash"] != second["data_hash"] or first["outer_hash"] != second["outer_hash"]:
            raise ValueError("Paired layer data or final outer fits differ: " + str(key))
        if abs(float(first["truth"])-float(second["truth"])) > 1e-12:
            raise ValueError("Paired layer truths differ")
        errors = [float(row["estimate"])-float(row["truth"]) for row in (first, second)]
        paired[key[:-1]].append(dict(errors=errors, variances=[float(row["variance"]) for row in (first,second)]))
        for layers, row in ((2, first), (3, second)):
            paired_groups[key[:-1]+(layers,)].append(row)
            cell, recipe, mode, quantity, seed = key
            if quantity == "tate":
                components = {arm: pair_records[(cell,recipe,mode,arm,seed)][layers] for arm in ("mu1","mu0")}
                paired_covariance_records[(cell,recipe,mode,layers)].append(components)
    paired_summaries = []
    for key, count in sorted(planned.items()):
        cell, recipe, mode, quantity, layers = key
        paired_summaries.append(dict(zip(CELL, cell), recipe=recipe, aggregation_mode=mode, quantity=quantity,
            crossfit_layers=layers, n_complete_in_layer=len(groups[key]), **metrics(paired_groups[key],count)))
    differences = []
    for key, values in sorted(paired.items()):
        cell, recipe, mode, quantity = key
        first = [value["errors"][0] for value in values]
        second = [value["errors"][1] for value in values]
        mse = [x*x-y*y for x,y in zip(first,second)]
        variances = paired_variance_comparison(first,second,[value["variances"][0] for value in values],
                                               [value["variances"][1] for value in values])
        differences.append(dict(zip(CELL, cell), recipe=recipe, aggregation_mode=mode, quantity=quantity,
            n_pairs=len(values), mean_estimate_difference_two_minus_three=statistics.mean(x-y for x,y in zip(first,second)),
            mse_difference_two_minus_three=statistics.mean(mse),
            mse_difference_mcse=statistics.stdev(mse)/math.sqrt(len(mse)) if len(mse)>1 else None, **variances))
    weight_groups = defaultdict(list)
    for key, values in weight_records.items():
        weight_groups[key[:-1]].append((statistics.mean(values), sum(abs(value)<=1e-10 for value in values)/len(values)))
    weights = []
    for (cell, recipe, mode, coordinate, layers), values in sorted(weight_groups.items()):
        weights.append(dict(zip(CELL, cell), recipe=recipe, aggregation_mode=mode, coordinate=coordinate,
            crossfit_layers=layers, n_repeats=len(values), mean_weight=statistics.mean(value[0] for value in values),
            sd_of_repeat_mean_weight=statistics.stdev(value[0] for value in values) if len(values)>1 else None,
            mean_fraction_outer_weights_near_zero=statistics.mean(value[1] for value in values)))
    write_union_csv(output/"metrics.csv", summaries)
    write_union_csv(output/"paired_metrics.csv", paired_summaries)
    write_union_csv(output/"paired_differences.csv", differences)
    write_union_csv(output/"arm_covariance.csv", covariance_summary(covariance_records))
    write_union_csv(output/"paired_arm_covariance.csv", covariance_summary(paired_covariance_records))
    write_union_csv(output/"paired_variance_components.csv", variance_component_summary(paired_covariance_records))
    write_union_csv(output/"weights.csv", weights)
    write_union_csv(output/"task_status.csv", statuses)
    write_union_csv(output/"runtime.csv", [dict(zip(CELL, cell), crossfit_layers=layers,
        n_success=len(values), median_seconds=statistics.median(values)) for (cell,layers),values in sorted(runtime.items())])
    text = ["# Two/three-layer comparison", "",
            "Committed repeats: {}/{}.".format(sum(row["committed"] for row in statuses),len(statuses)), "",
            "metrics.csv includes both arm means and TATE, empirical and estimated variances, fixed-weight inference and coverage intervals.",
            "paired_metrics.csv restricts those metrics to the identical completed seed pairs; use it for interim layer comparisons.",
            "paired_arm_covariance.csv uses that same paired subset for both arm covariance estimates.",
            "paired_variance_components.csv decomposes TATE variance into both arm variances and minus twice their covariance, with delete-one jackknife uncertainty for reported-minus-empirical discrepancies.",
            "paired_differences.csv compares only identical data and outer-fit fingerprints; unavailable pairs remain visible in task_status.csv.",
            "The variance discrepancy compares mean reported variance differences with empirical variance differences; its jackknife SE uses paired repeats.",
            "Eta sensitivity holds nuisance fits and active sets fixed. It is not full-refit bootstrap inference.",
            "Runtime is descriptive: machines, checkpoint reuse and restarts can differ. Near-zero weights use absolute tolerance 1e-10.",
            "Small cells have substantial Monte Carlo uncertainty; no seeds are filtered by coverage or error."]
    (output/"report.md").write_text("\n".join(text)+"\n")
    print(text[2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    review(args.plan_root,args.output)

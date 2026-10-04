"""Prespecified paired residual-inference experiments; artifacts only on scratch."""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
import numpy as np
from scipy.stats import chi2, norm
from inference import (TargetRegion, count_interval_set, dispersion_intervals,
                       fourier_interval_set, fourier_summary, intersect_sets)

CORRELATION = -.25
MAX_SOURCES = 2048


def settings():
    result = []
    def add(panel, counts, patterns, floors, noise="gaussian", **extra):
        for count in counts:
            for pattern in patterns:
                for floor in floors:
                    row = dict(panel=panel, source_count=count, pattern=pattern,
                               shared_scale=floor, noise=noise, valid_fraction=.75,
                               heterogeneous=False, n_per_site=1000, frequency_scale=2.)
                    row.update(extra)
                    row["source_evaluation_size"] = row["n_per_site"] - row.get("pilot_size", 0)
                    row["target_evaluation_size"] = row["n_per_site"]
                    row["cell_id"] = len(result) + 1
                    result.append(row)
    add("main", [8, 32, 128, 512, 2048], ["all_valid", "weak_below", "weak_boundary", "strong"], [0., .5])
    add("bias_shape", [32, 128, 512], ["weak_disjoint", "arm_cancel", "one_extreme"], [.5])
    add("heterogeneous", [32, 128, 512], ["all_valid", "weak_boundary", "strong"], [.5], heterogeneous=True)
    for fraction in [.5, .25]:
        add("valid_count", [128, 512], ["all_valid", "weak_boundary"], [.5], valid_fraction=fraction)
    add("variance_pilot", [32, 512, 2048], ["all_valid", "weak_boundary"], [.5], pilot_size=250)
    add("binary_scores", [32, 512, 2048], ["all_valid", "weak_boundary", "strong"], [0., .3], noise="binary")
    return result


def source_bias(setting, standard_error):
    count = setting["source_count"]
    pattern = setting["pattern"]
    if pattern == "ambiguous_centers":
        bias = np.zeros((count, 2))
        valid_per_four = int(round(4 * setting["actual_valid_fraction"]))
        bias[np.arange(count) % 4 >= valid_per_four, 0] = standard_error
        return bias
    scales = {"all_valid": 0., "weak_below": 2 * count**(-1/3),
              "weak_boundary": 2 * count**(-1/4), "strong": 4.,
              "weak_disjoint": 2 * count**(-1/4), "arm_cancel": 4., "one_extreme": 12.}
    bias = np.zeros((count, 2))
    if pattern == "one_extreme":
        bias[0, 0] = scales[pattern] * standard_error
    else:
        bias[np.arange(count) % 4 == 0, 0] = scales[pattern] * standard_error
        if pattern == "weak_disjoint":
            bias[np.arange(count) % 4 == 1, 1] = -scales[pattern] * standard_error
        elif pattern == "arm_cancel":
            bias[:, 1] = bias[:, 0]
    return bias


def gaussian_data(setting, repeats):
    count = setting["source_count"]
    evaluation_size = setting.get("source_evaluation_size", setting["n_per_site"])
    sigma = math.sqrt(.5 / evaluation_size)
    target_sigma = math.sqrt(.5 / setting["n_per_site"])
    # Common random numbers across Gaussian K and bias settings. Site identities
    # are nested, with every fourth source invalid in the treated arm.
    rng = np.random.RandomState(940117 + setting.get("seed_offset", 0))
    independent = rng.normal(size=(repeats, MAX_SOURCES, 2))[:, :count]
    independent[:, :, 1] = CORRELATION * independent[:, :, 0] + math.sqrt(1 - CORRELATION**2) * independent[:, :, 1]
    factors = np.ones((count, 2))
    if setting["heterogeneous"]:
        factors[:, 0] = .7 + .6 * (np.arange(count) % 7) / 6
        factors[:, 1] = 1.3 - .6 * (np.arange(count) % 7) / 6
    standard_errors = sigma * factors
    variances = standard_errors**2
    bias = source_bias(setting, sigma)
    residual_truth = np.array([.02, -.01])
    source = residual_truth + bias + independent * standard_errors
    standard_deviations = target_sigma * np.array([setting["shared_scale"], 1., 1.])
    target_correlation = np.array([[1, .2, -.15], [.2, 1, CORRELATION], [-.15, CORRELATION, 1.]])
    target_covariance = target_correlation * np.outer(standard_deviations, standard_deviations)
    target_rng = np.random.RandomState(713905 + setting.get("seed_offset", 0))
    target_mean = np.array([.1, .02, -.01])
    target = target_mean + (target_rng.normal(size=(repeats, 3)) @ np.linalg.cholesky(target_correlation).T) * standard_deviations
    data = dict(source=source, variance=variances, bias=bias, residual_truth=residual_truth,
                target=target, target_covariance=target_covariance, truth=.13,
                source_cross_covariance=CORRELATION * np.prod(standard_errors, axis=1))
    if setting["panel"] == "variance_pilot":
        degrees = setting["pilot_size"] - 1
        rng = np.random.RandomState(608933)
        first = np.sqrt(rng.chisquare(degrees, size=(repeats, count)))
        off_diagonal = rng.normal(size=(repeats, count))
        second = np.sqrt(rng.chisquare(degrees - 1, size=(repeats, count)))
        estimates = np.empty_like(source)
        estimates[:, :, 0] = variances[:, 0] * first**2 / degrees
        estimates[:, :, 1] = variances[:, 1] * ((CORRELATION * first + math.sqrt(1 - CORRELATION**2) * off_diagonal)**2 +
            (1 - CORRELATION**2) * second**2) / degrees
        data["pilot_variance"] = estimates
    return data


def binary_data(setting, repeats):
    """Exact multinomial patient simulation: X=+/-1, A randomized, binary Y.

    m1=.55+cX/2, m0=.45-cX/2; target truth=.10 and both true residuals=0.
    Oracle merged weights equal 2. Each source has exactly 1000 patients.
    """
    count, heterogeneity = setting["source_count"], setting["shared_scale"]
    size = setting.get("source_evaluation_size", setting["n_per_site"])
    target_size = setting.get("target_evaluation_size", setting["n_per_site"])
    pilot_size = setting.get("pilot_size", 0)
    base_variance = .495 - heterogeneity**2 / 2
    bias = source_bias(setting, math.sqrt(base_variance / size))
    states = np.array([(x, arm, outcome) for x in [-1, 1] for arm in [0, 1] for outcome in [0, 1]])
    prediction = np.column_stack((.55 + heterogeneity * states[:, 0] / 2,
                                  .45 - heterogeneity * states[:, 0] / 2))
    residual = np.column_stack((2 * (states[:, 1] == 1) * (states[:, 2] - prediction[:, 0]),
                                2 * (states[:, 1] == 0) * (states[:, 2] - prediction[:, 1])))
    def state_probabilities(shift):
        probability = np.where(states[:, 1] == 1, prediction[:, 0] + shift[0], prediction[:, 1] + shift[1])
        if np.any((probability <= 0) | (probability >= 1)):
            raise ValueError("Binary scenario violates outcome probability range")
        return .25 * np.where(states[:, 2] == 1, probability, 1 - probability)
    target_scores = np.column_stack((prediction[:, 0] - prediction[:, 1], residual))
    target_probability = state_probabilities(np.zeros(2))
    target_mean = target_probability @ target_scores
    target_centered = target_scores - target_mean
    target_covariance = target_centered.T @ (target_probability[:, None] * target_centered) / target_size
    target_rng = np.random.RandomState(413905 + setting.get("seed_offset", 0))
    target_counts = target_rng.multinomial(target_size, target_probability, size=repeats)
    target = target_counts @ target_scores / target_size
    target_second = np.einsum("bs,si,sj->bij", target_counts, target_scores, target_scores)
    target_sample_covariance = (target_second - target_size * np.einsum("bi,bj->bij", target, target)) / (target_size * (target_size - 1))
    source = np.empty((repeats, count, 2))
    sample_variance = np.empty_like(source)
    variances = np.empty((count, 2))
    cross_covariance = np.empty(count)
    # Exact patient-state counts avoid materializing all individual records.
    rng = np.random.RandomState(510000 + count + setting.get("seed_offset", 0))
    pilot_rng = np.random.RandomState(610000 + count + setting.get("seed_offset", 0))
    pilot_variance = np.empty_like(source) if pilot_size else None
    for shift in np.unique(bias, axis=0):
        indices = np.flatnonzero(np.all(bias == shift, axis=1))
        probability = state_probabilities(shift)
        counts = rng.multinomial(size, probability, size=(repeats, len(indices)))
        means = counts @ residual / size
        second = counts @ (residual**2)
        source[:, indices] = means
        sample_variance[:, indices] = (second - size * means**2) / (size * (size - 1))
        if pilot_size:
            pilot_counts = pilot_rng.multinomial(pilot_size, probability, size=(repeats, len(indices)))
            pilot_means = pilot_counts @ residual / pilot_size
            pilot_second = pilot_counts @ (residual**2)
            pilot_variance[:, indices] = (pilot_second - pilot_size * pilot_means**2) / ((pilot_size - 1) * size)
        expectation = probability @ residual
        centered = residual - expectation
        covariance = centered.T @ (probability[:, None] * centered) / size
        variances[indices] = np.diag(covariance)
        cross_covariance[indices] = covariance[0, 1]
    return dict(source=source, variance=variances, bias=bias, residual_truth=np.zeros(2),
                target=target, target_covariance=target_covariance, truth=.1,
                target_sample_covariance=target_sample_covariance, target_support=target_scores,
                target_sample_size=target_size, pilot_variance=pilot_variance,
                source_cross_covariance=cross_covariance, sample_variance=sample_variance,
                score_absolute_bounds=np.max(abs(residual), axis=0),
                score_ranges=np.ptp(residual, axis=0), sample_sizes=np.full(count, size))


def compute_cell(setting, repeats, output_directory):
    started = time.time()
    if setting.get("hybrid", False) and setting["noise"] != "gaussian":
        raise NotImplementedError("Combined-certificate experiment currently requires Gaussian oracle scores")
    count = setting["source_count"]
    valid_minimum = int(math.ceil(setting["valid_fraction"] * count))
    data = gaussian_data(setting, repeats) if setting["noise"] == "gaussian" else binary_data(setting, repeats)
    source, variances, target = data["source"], data["variance"], data["target"]
    target_covariance = data["target_covariance"]
    truth = data["truth"]
    target_region = TargetRegion(target_covariance)
    moment_region = TargetRegion(target_covariance, calibration="moment")
    contrast = np.array([1., 1., -1.])
    target_point = target @ contrast
    target_radius = norm.isf(.025) * math.sqrt(target_region.effect_variance)
    methods = {}
    def store(name, interval, point=None, fallback=None, accepted=None, branches=None):
        methods[name] = dict(interval=np.asarray(interval), point=np.asarray(point) if point is not None else np.mean(interval, axis=1),
            fallback=np.zeros(repeats, dtype=bool) if fallback is None else fallback,
            source_accept=np.full((repeats, 2), np.nan) if accepted is None else accepted,
            branches=np.full((repeats, 2), np.nan) if branches is None else branches)
    store("target_nominal", np.column_stack((target_point - target_radius, target_point + target_radius)), target_point)
    joint_radius = math.sqrt(target_region.critical_squared * target_region.effect_variance)
    store("target_joint", np.column_stack((target_point - joint_radius, target_point + joint_radius)), target_point)
    source_point = target[:, 0] + source[:, :, 0].mean(axis=1) - source[:, :, 1].mean(axis=1)
    source_variance = target_covariance[0, 0] + (np.sum(variances) - 2 * np.sum(data["source_cross_covariance"])) / count**2
    pool_radius = norm.isf(.025) * math.sqrt(source_variance)
    store("pooled_naive", np.column_stack((source_point - pool_radius, source_point + pool_radius)), source_point)
    valid = data["bias"] == 0
    valid_counts = valid.sum(axis=0)
    oracle_point = target[:, 0] + source[:, valid[:, 0], 0].mean(axis=1) - source[:, valid[:, 1], 1].mean(axis=1)
    oracle_variance = target_covariance[0, 0] + sum(np.sum(variances[valid[:, arm], arm]) / valid_counts[arm]**2 for arm in [0, 1]) - \
        2 * np.sum(data["source_cross_covariance"][np.all(valid, axis=1)]) / np.prod(valid_counts)
    oracle_radius = norm.isf(.025) * math.sqrt(oracle_variance)
    store("oracle_valid", np.column_stack((oracle_point - oracle_radius, oracle_point + oracle_radius)), oracle_point)

    summaries = {}
    for calibration in ["hoeffding", "bernstein"]:
        summaries["fourier_" + calibration] = [fourier_summary(source[:, :, arm], variances[:, arm], .0125,
            setting["frequency_scale"], calibration=calibration) for arm in [0, 1]]
    if setting.get("hybrid", False):
        # Split each source-arm error allowance before observing confirmation
        # data. The two certificates need not be independent.
        summaries["hybrid"] = [fourier_summary(source[:, :, arm], variances[:, arm], .00625,
            setting["frequency_scale"]) for arm in [0, 1]]
        hybrid_dispersion = [dispersion_intervals(source[:, :, arm], variances[:, arm],
            valid_minimum, .00625) for arm in [0, 1]]
    if setting["panel"] == "variance_pilot":
        pilot = data["pilot_variance"]
        summaries["fourier_plugin"] = [fourier_summary(source[:, :, arm], pilot[:, :, arm], .0125,
            setting["frequency_scale"]) for arm in [0, 1]]
        degrees = setting["pilot_size"] - 1
        tail = .005 / (4 * count)
        lower = pilot * degrees / chi2.isf(tail, degrees)
        upper = pilot * degrees / chi2.ppf(tail, degrees)
        summaries["fourier_variance_certified"] = [fourier_summary(source[:, :, arm], lower[:, :, arm], .01,
            setting["frequency_scale"], variance_upper=upper[:, :, arm]) for arm in [0, 1]]
    if setting["noise"] == "binary":
        summaries["fourier_plugin"] = [fourier_summary(source[:, :, arm],
            np.maximum(data["sample_variance"][:, :, arm], 1e-14), .0125, setting["frequency_scale"]) for arm in [0, 1]]
        summaries["fourier_bounded"] = [fourier_summary(source[:, :, arm], variances[:, arm], .0125,
            setting["frequency_scale"], score_ranges=data["score_ranges"][arm],
            sample_sizes=data["sample_sizes"]) for arm in [0, 1]]
    for method, arm_summaries in summaries.items():
        region = moment_region if method == "fourier_bounded" else target_region
        intervals = np.empty((repeats, 2))
        fallback = np.zeros(repeats, dtype=bool)
        accepted = np.empty((repeats, 2), dtype=bool)
        branches = np.zeros((repeats, 2), dtype=int)
        for repetition in range(repeats):
            domains = region.domains(target[repetition])
            source_sets = []
            for arm in [0, 1]:
                transform, frequencies, slack = (item[repetition] for item in arm_summaries[arm])
                source_sets.append(fourier_interval_set(transform, frequencies, slack, valid_minimum / count, domains[arm]))
                accepted[repetition, arm] = np.all(abs(transform - valid_minimum / count *
                    np.exp(1j * frequencies * data["residual_truth"][arm])) <= 1 - valid_minimum / count + slack)
                if method == "hybrid":
                    dispersion_interval = hybrid_dispersion[arm][repetition]
                    source_sets[-1] = intersect_sets(source_sets[-1], [tuple(dispersion_interval)])
                    accepted[repetition, arm] &= dispersion_interval[0] <= data["residual_truth"][arm] <= dispersion_interval[1]
                branches[repetition, arm] = len(source_sets[arm])
            projected = region.project(target[repetition], source_sets)
            fallback[repetition] = projected is None
            intervals[repetition] = methods["target_nominal"]["interval"][repetition] if projected is None else projected
        store(method, intervals, fallback=fallback, accepted=accepted, branches=branches)

    bounded = setting["noise"] == "binary"
    dispersion = [dispersion_intervals(source[:, :, arm], variances[:, arm], valid_minimum, .0125,
        sample_variances=data["sample_variance"][:, :, arm] if bounded else None,
        sample_sizes=data.get("sample_sizes"),
        score_ranges=data["score_ranges"][arm] if bounded else None) for arm in [0, 1]]
    for method in ["dispersion", "count_search"]:
        region = moment_region if bounded and method == "dispersion" else target_region
        intervals = np.empty((repeats, 2))
        fallback = np.zeros(repeats, dtype=bool)
        accepted = np.empty((repeats, 2), dtype=bool)
        branches = np.zeros((repeats, 2), dtype=int)
        for repetition in range(repeats):
            domains = region.domains(target[repetition])
            if method == "dispersion":
                source_sets = [[tuple(dispersion[arm][repetition])] for arm in [0, 1]]
                accepted[repetition] = [dispersion[arm][repetition, 0] <= data["residual_truth"][arm] <=
                                        dispersion[arm][repetition, 1] for arm in [0, 1]]
            else:
                source_sets = [count_interval_set(source[repetition, :, arm], variances[:, arm], valid_minimum, .0125,
                               domains[arm]) for arm in [0, 1]]
                accepted[repetition] = [any(lower <= data["residual_truth"][arm] <= upper for lower, upper in source_sets[arm])
                                       for arm in [0, 1]]
            branches[repetition] = [len(item) for item in source_sets]
            projected = region.project(target[repetition], source_sets)
            fallback[repetition] = projected is None
            intervals[repetition] = methods["target_nominal"]["interval"][repetition] if projected is None else projected
        store(method, intervals, fallback=fallback, accepted=accepted, branches=branches)

    names = list(methods)
    lower = np.column_stack([methods[name]["interval"][:, 0] for name in names])
    upper = np.column_stack([methods[name]["interval"][:, 1] for name in names])
    points = np.column_stack([methods[name]["point"] for name in names])
    fallbacks = np.column_stack([methods[name]["fallback"] for name in names])
    if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(lower > upper):
        raise ArithmeticError("Nonfinite or reversed confidence interval")
    metrics = []
    for column, name in enumerate(names):
        covered = (lower[:, column] <= truth) & (truth <= upper[:, column])
        lengths = upper[:, column] - lower[:, column]
        source_accept = methods[name]["source_accept"]
        rows = dict(setting, method=name, repeats=repeats, truth=truth,
                    coverage=float(covered.mean()), coverage_mcse=float(np.sqrt(covered.mean() * (1 - covered.mean()) / repeats)),
                    mean_length=float(lengths.mean()), length_mcse=float(lengths.std(ddof=1) / math.sqrt(repeats)),
                    length_ratio_target=float(lengths.mean() / (2 * target_radius)),
                    bias=float(np.mean(points[:, column] - truth)), rmse=float(np.sqrt(np.mean((points[:, column] - truth)**2))),
                    fallback_rate=float(fallbacks[:, column].mean()),
                    source_mu1_accept=float(np.mean(source_accept[:, 0])), source_mu0_accept=float(np.mean(source_accept[:, 1])),
                    max_branches=float(np.nanmax(methods[name]["branches"])) if np.any(np.isfinite(methods[name]["branches"])) else float("nan"))
        metrics.append(rows)
    path = Path(output_directory) / "cells" / ("cell_{:03d}".format(setting["cell_id"]))
    np.savez_compressed(str(path) + ".npz", methods=np.array(names), lower=lower, upper=upper,
                        point=points, fallback=fallbacks, source_bias=data["bias"], source_variance=variances,
                        source_accept=np.stack([methods[name]["source_accept"] for name in names], axis=1),
                        source_branches=np.stack([methods[name]["branches"] for name in names], axis=1),
                        target_covariance=target_covariance)
    (Path(str(path) + ".json")).write_text(json.dumps(dict(setting=setting, metrics=metrics,
        runtime_seconds=time.time() - started), indent=2) + "\n")
    return metrics, time.time() - started


def run_task(arguments):
    return compute_cell(*arguments)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--repeats", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=2)
    profiles = parser.add_mutually_exclusive_group()
    profiles.add_argument("--smoke", action="store_true")
    profiles.add_argument("--confirmation", action="store_true")
    profiles.add_argument("--boundary", action="store_true")
    profiles.add_argument("--variance-budget", action="store_true")
    arguments = parser.parse_args()
    if arguments.repeats < 1 or not 1 <= arguments.workers <= 4:
        raise ValueError("Positive repeat count and one to four allocated workers required")
    output = Path(arguments.output)
    if not str(output).startswith("/scratch.global/zhan9381/FACE-HD/") or output.exists():
        raise ValueError("Use a new FACE-HD scratch output directory")
    output.mkdir(parents=True)
    (output / "cells").mkdir()
    plan = settings()
    if arguments.variance_budget:
        plan = [dict(row, source_evaluation_size=row["n_per_site"] - row["pilot_size"])
                for row in plan if row["panel"] == "variance_pilot"]
    if arguments.boundary:
        template = plan[0]
        plan = [dict(template, panel="minority_boundary", source_count=count,
                     pattern="ambiguous_centers", valid_fraction=fraction,
                     actual_valid_fraction=fraction, hybrid=True, seed_offset=3000000,
                     cell_id=index + 1)
                for index, (count, fraction) in enumerate((count, fraction)
                    for count in [32, 512, 2048] for fraction in [.5, .25])]
    if arguments.confirmation:
        plan = [dict(row, panel="confirmation", hybrid=True, seed_offset=2000000)
                for row in plan if row["panel"] == "main" and row["source_count"] >= 128
                and row["pattern"] in ["all_valid", "weak_boundary", "strong"]]
    if arguments.smoke:
        plan = [plan[0], next(row for row in plan if row["panel"] == "heterogeneous"),
                next(row for row in plan if row["panel"] == "variance_pilot"),
                next(row for row in plan if row["panel"] == "binary_scores")]
    source_directory = Path(__file__).resolve().parent
    hashes = {file.name: hashlib.sha256(file.read_bytes()).hexdigest() for file in source_directory.glob("*.py")}
    configuration = dict(settings=plan, repeats=arguments.repeats, workers=arguments.workers,
        code_sha256=hashes, method_scope="Oracle residual-score experiments, not fitted high-dimensional RoCE.",
        alpha_target=.025, alpha_source=.0125, certified_variance_alpha=.005,
        certified_variance_source_alpha=.01, gaussian_pairing="Common patient-score draws across K/bias settings; methods paired in every cell.",
        point_rule="Interval midpoint for compatibility methods; ordinary raw estimates for baselines.")
    (output / "configuration.json").write_text(json.dumps(configuration, indent=2) + "\n")
    tasks = [(row, arguments.repeats, str(output)) for row in plan]
    all_metrics = []
    with multiprocessing.Pool(arguments.workers) as pool:
        for metrics, runtime in pool.imap_unordered(run_task, tasks):
            all_metrics.extend(metrics)
            print("Completed cell {} ({}, K={}) in {:.1f}s".format(metrics[0]["cell_id"],
                  metrics[0]["panel"], metrics[0]["source_count"], runtime), flush=True)
    fields = sorted(set(key for row in all_metrics for key in row))
    with (output / "metrics.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(all_metrics, key=lambda row: (row["cell_id"], row["method"])))
    (output / "COMPLETE").write_text("All prespecified cells completed; every repetition retained.\n")


if __name__ == "__main__":
    main()

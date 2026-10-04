"""Protected fallback, explicit budgets and directional target calibration.

These research components extend inference.py without changing frozen studies.
Every confidence allowance is explicit; oracle and empirical inputs are distinct.
"""
import math
import numpy as np
from scipy.stats import norm
from concentrated_bounds import concentrated_dispersion_interval
from inference import (TargetRegion, fourier_summary, fourier_interval_set,
                       dispersion_intervals, intersect_sets)

DIRECTIONS = np.array([[1., 0, 0], [0, 1., 0], [0, 0, 1.], [1., 1., -1.]])


def directional_radii(covariance, score_ranges, sample_size, alpha, calibration,
                       coordinate_bias=None):
    """Four simultaneous directions, including the complete target TATE.

    empirical: Maurer--Pontil Theorem 4, two-sided, with unbiased mean covariance.
    bernstein: known/certified mean variance and deterministic patient ranges.
    normal: diagnostic marginal Gaussian calibration, not a fitted-score theorem.
    """
    if calibration not in ["normal", "bernstein", "empirical"] or not 0 < alpha < 1:
        raise ValueError("Unknown calibration or invalid error allowance")
    covariance = np.asarray(covariance)
    variances = np.maximum(0, np.einsum("di,...ij,dj->...d", DIRECTIONS, covariance, DIRECTIONS))
    if calibration == "normal":
        radius = norm.isf(alpha / 8) * np.sqrt(variances)
    else:
        if score_ranges is None or sample_size < 2:
            raise ValueError("Bounded calibrations require patient ranges and sample size")
        ranges = np.asarray(score_ranges)
        logarithm = math.log((16 if calibration == "empirical" else 8) / alpha)
        if calibration == "empirical":
            radius = np.sqrt(2 * variances * logarithm) + 7 * ranges * logarithm / (3 * (sample_size - 1))
        else:
            term = ranges * logarithm / (3 * sample_size)
            radius = term + np.sqrt(term**2 + 2 * variances * logarithm)
        radius = np.where(ranges == 0, 0, radius)
    if coordinate_bias is not None:
        radius = radius + np.asarray(coordinate_bias) @ abs(DIRECTIONS).T
    return radius


class DirectionalTargetRegion:
    """Three coordinate bands plus a covariance-aware TATE band."""
    def __init__(self, radii):
        self.radii = np.asarray(radii)
        if self.radii.shape != (4,) or np.any(self.radii < 0) or not np.all(np.isfinite(self.radii)):
            raise ValueError("Four finite nonnegative directional radii required")

    def domains(self, center):
        return [(center[index] - self.radii[index], center[index] + self.radii[index]) for index in [1, 2]]

    def project(self, center, residual_sets):
        domains = self.domains(center)
        first = intersect_sets(residual_sets[0], [domains[0]])
        second = intersect_sets(residual_sets[1], [domains[1]])
        effect = float(DIRECTIONS[3] @ center)
        target_band = (effect - self.radii[3], effect + self.radii[3])
        intervals = []
        for treated in first:
            for control in second:
                lower = max(center[0] - self.radii[0] + treated[0] - control[1], target_band[0])
                upper = min(center[0] + self.radii[0] + treated[1] - control[0], target_band[1])
                if lower <= upper + 1e-12:
                    intervals.append((lower - 1e-11, max(lower, upper) + 1e-11))
        if not intervals:
            return None
        return min(item[0] for item in intervals), max(item[1] for item in intervals)


def variance_certificate(pilot_mean_variance, ranges, pilot_size, evaluation_size, alpha):
    """Simultaneous bounded-score standard-deviation bounds, MP Theorem 10.

    Last two axes index sites and arms; pilot_mean_variance is S_pilot^2/n_eval.
    Two tails for each of 2K variances share alpha. Deterministic range variance
    limits cap the upper endpoint. Pilot observations are independent of fits
    and final evaluation observations.
    """
    estimate = np.asarray(pilot_mean_variance)
    count = estimate.shape[-2]
    if pilot_size < 2 or evaluation_size < 2 or not 0 < alpha < 1:
        raise ValueError("Invalid variance-certificate sample sizes or allowance")
    radius = np.asarray(ranges) * math.sqrt(2 * math.log(4 * count / alpha) /
                                         ((pilot_size - 1) * evaluation_size))
    deviation = np.sqrt(np.maximum(estimate, 0))
    lower = np.maximum(0, deviation - radius)**2
    upper = np.minimum((deviation + radius)**2, np.asarray(ranges)**2 / (4 * evaluation_size))
    # The sample SD can exceed the population range bound by n/(n-1). Capping
    # the lower endpoint also avoids impossible numerical certificate ordering.
    lower = np.minimum(lower, upper)
    return lower, upper


class SourceCompatibilityModule:
    """One configurable source procedure for the paired component ablations."""
    def __init__(self, data, valid_minimum, alpha_arm, variance_lower=None,
                 variance_upper=None, bounded=False, source_mode="hybrid",
                 dispersion_calibration="cantelli", source_domain=(-1., 1.)):
        if source_mode not in ["hybrid", "fourier", "dispersion"]:
            raise ValueError("Unknown source compatibility mode")
        if dispersion_calibration not in ["cantelli", "bounded_exponential"]:
            raise ValueError("Unknown dispersion calibration")
        if dispersion_calibration == "bounded_exponential" and not bounded:
            raise ValueError("Exponential dispersion requires bounded patient scores")
        if source_mode == "fourier" and not bounded:
            raise ValueError("Fourier-only prototype requires a declared bounded residual domain")
        self.mode = source_mode
        self.domain = source_domain
        self.count = data["source"].shape[1]
        counts = np.broadcast_to(valid_minimum, (2,))
        if not np.all(np.isfinite(counts)) or np.any(counts != np.floor(counts)):
            raise ValueError("Valid-source counts must be finite integers")
        self.valid_minimum = counts.astype(int)
        self.fourier = []
        self.moments = []
        lower = data["variance"] if variance_lower is None else variance_lower
        upper = lower if variance_upper is None else variance_upper
        alpha_component = alpha_arm / 2 if source_mode == "hybrid" else alpha_arm
        if np.any(self.valid_minimum < 1) or np.any(self.valid_minimum > self.count):
            raise ValueError("At least one guaranteed source in each arm required")
        for arm in [0, 1]:
            score_range = np.asarray(data["score_ranges"])[..., arm] if bounded else None
            sizes = data.get("sample_sizes") if bounded or variance_upper is not None else None
            sample_variance = data["sample_variance"][:, :, arm] if bounded or variance_upper is not None else None
            if source_mode != "dispersion":
                self.fourier.append(fourier_summary(data["source"][:, :, arm], np.asarray(lower)[..., arm],
                    alpha_component, variance_upper=np.asarray(upper)[..., arm] if variance_upper is not None else None,
                    score_ranges=score_range, sample_sizes=sizes))
            if source_mode == "fourier":
                self.moments.append(np.broadcast_to(source_domain, (len(data["source"]), 2)))
            elif dispersion_calibration == "bounded_exponential":
                self.moments.append(concentrated_dispersion_interval(data["source"][:, :, arm],
                    sample_variance, sizes, score_range, np.asarray(data["score_absolute_bounds"])[..., arm],
                    self.valid_minimum[arm], alpha_component))
            else:
                self.moments.append(dispersion_intervals(data["source"][:, :, arm], np.asarray(upper)[..., arm],
                    self.valid_minimum[arm], alpha_component, sample_variances=sample_variance,
                    sample_sizes=sizes, score_ranges=score_range,
                    variance_status="certified" if variance_upper is not None or bounded else "known_gaussian"))

    def intervals(self, repetition, budgets, squared_budgets=None):
        residual_sets = []
        moment_intervals = []
        for arm in [0, 1]:
            base = self.moments[arm][repetition]
            interval = self.domain if self.mode == "fourier" else (base[0] - budgets[arm], base[1] + budgets[arm])
            moment_intervals.append(interval)
            if self.mode == "dispersion":
                residual_sets.append([interval])
                continue
            transform, frequencies, allowance = (item[repetition] for item in self.fourier[arm])
            fraction = self.valid_minimum[arm] / self.count
            drift_allowance = fraction * abs(frequencies) * budgets[arm]
            if squared_budgets is not None:
                # budgets then control the signed mean over an unknown valid
                # g-subset; squared_budgets bound its squared drifts divided by K.
                drift_allowance += frequencies**2 * squared_budgets[arm] / 2
            allowance = allowance + np.minimum(2 * fraction, drift_allowance)
            residual_sets.append(fourier_interval_set(transform, frequencies, allowance, fraction, interval))
        return residual_sets, moment_intervals


def protected_projection(center, target_region, source_sets, moment_intervals,
                         prediction_radius, target_reference, fallback_rule="shorter"):
    """Return interval, source-compatible point, fallback code and length envelope.

    On the simultaneous good event the target/source intersection is nonempty.
    A fallback cannot remove that event. The shorter rule additionally ensures
    every realized interval is no wider than the standalone moment envelope.
    Codes: 0=nonempty projection, 1=target fallback, 2=source fallback.
    """
    if fallback_rule not in ["target", "shorter"]:
        raise ValueError("Unknown fallback rule")
    moment_band = (center[0] - prediction_radius + moment_intervals[0][0] - moment_intervals[1][1],
                   center[0] + prediction_radius + moment_intervals[0][1] - moment_intervals[1][0])
    if all(source_sets):
        hulls = [(intervals[0][0], intervals[-1][1]) for intervals in source_sets]
    else:
        hulls = moment_intervals
    source_band = (center[0] - prediction_radius + hulls[0][0] - hulls[1][1],
                   center[0] + prediction_radius + hulls[0][1] - hulls[1][0])
    source_hint = center[0] + np.mean(hulls[0]) - np.mean(hulls[1])
    interval = target_region.project(center, source_sets)
    code = 0
    if interval is None:
        interval, code = target_reference, 1
        if fallback_rule == "shorter" and source_band[1] - source_band[0] < interval[1] - interval[0]:
            interval, code = source_band, 2
    width = interval[1] - interval[0]
    envelope_width = moment_band[1] - moment_band[0]
    if fallback_rule == "shorter" and width > envelope_width + 1e-8:
        raise ArithmeticError("Protected fallback exceeded its pointwise length envelope")
    return interval, float(np.clip(source_hint, *interval)), code, envelope_width


def infer_repair(data, target_alpha=.025, fallback_rule="shorter", target_calibration="gaussian",
                 variance_lower=None, variance_upper=None, variance_failure=0., bias_failure=0.,
                 source_bias_budget=None, target_bias_budget=None, valid_minimum=None, bounded=False,
                 source_mode="hybrid", dispersion_calibration="cantelli",
                 source_squared_bias_budget=None, source_domain=(-1., 1.)):
    if data.get("source_evaluation_reuses_design", False) and (source_mode != "dispersion" or dispersion_calibration != "bounded_exponential"):
        raise ValueError("Reused source designs are supported only by direct bounded dispersion")
    repeats, count = data["source"].shape[:2]
    if target_calibration not in ["gaussian", "moment", "normal", "bernstein", "empirical"]:
        raise ValueError("Unknown target region calibration")
    if valid_minimum is None:
        valid_minimum = math.ceil(.75 * count)
    if not 0 < target_alpha < .05 or not 0 <= variance_failure < .05 or not 0 <= bias_failure < .05:
        raise ValueError("Invalid target or certificate error allocation")
    alpha_arm = (.05 - target_alpha - variance_failure - bias_failure) / 2
    if alpha_arm <= 0:
        raise ValueError("Confidence allocations exceed .05")
    source_budget = np.zeros((repeats, 2)) if source_bias_budget is None else np.broadcast_to(source_bias_budget, (repeats, 2))
    target_budget = np.zeros((repeats, 3)) if target_bias_budget is None else np.broadcast_to(target_bias_budget, (repeats, 3))
    if any(not np.all(np.isfinite(values)) or np.any(values < 0) for values in [source_budget,target_budget]):
        raise ValueError("Bias budgets must be finite and nonnegative")
    squared_budget = None if source_squared_bias_budget is None else np.broadcast_to(source_squared_bias_budget, (repeats, 2))
    if squared_budget is not None and (np.any(squared_budget < 0) or not np.all(np.isfinite(squared_budget))):
        raise ValueError("Squared drift budgets must be finite and nonnegative")
    sources = SourceCompatibilityModule(data, valid_minimum, alpha_arm, variance_lower, variance_upper, bounded,
                                        source_mode, dispersion_calibration, source_domain)
    covariance = np.broadcast_to(data["target_covariance"], (repeats, 3, 3))
    mean = data["target"]
    effect = mean @ DIRECTIONS[3]
    target_variance = np.einsum("i,bij,j->b", DIRECTIONS[3], covariance, DIRECTIONS[3])
    reference_radius = norm.isf(.025) * np.sqrt(target_variance) + target_budget @ abs(DIRECTIONS[3])
    radii = None
    if target_calibration in ["normal", "bernstein", "empirical"]:
        chosen_covariance = data["target_sample_covariance"] if target_calibration in ["normal", "empirical"] else covariance
        if target_calibration == "normal":
            ranges = None
        elif "target_direction_ranges" in data:
            ranges = data["target_direction_ranges"]
        else:
            ranges = np.ptp(data["target_support"] @ DIRECTIONS.T, axis=0)
        radii = directional_radii(chosen_covariance, ranges, data["target_sample_size"], target_alpha,
                                  target_calibration, target_budget)
        radii = np.broadcast_to(radii, (repeats, 4))
        if target_calibration != "normal":
            # A valid, conservative target reference using the same target
            # directional calibration; the fallback does not require a new alpha.
            reference_radius = radii[:, 3]
    elif np.any(target_budget):
        raise ValueError("Biased target scores require a directional bias-aware region")
    intervals = np.empty((repeats, 2))
    points = np.empty(repeats)
    codes = np.empty(repeats, dtype=int)
    envelopes = np.empty(repeats)
    cached_region = None
    if radii is None and np.asarray(data["target_covariance"]).ndim == 2:
        cached_region = TargetRegion(data["target_covariance"], target_alpha,
                                     "moment" if target_calibration == "moment" else "gaussian")
    for repetition in range(repeats):
        if radii is not None:
            region = DirectionalTargetRegion(radii[repetition])
            prediction_radius = radii[repetition, 0]
        else:
            region = cached_region or TargetRegion(covariance[repetition], target_alpha,
                        "moment" if target_calibration == "moment" else "gaussian")
            prediction_radius = math.sqrt(region.critical_squared * covariance[repetition, 0, 0])
        residual_sets, moments = sources.intervals(repetition, source_budget[repetition],
            None if squared_budget is None else squared_budget[repetition])
        reference = (effect[repetition] - reference_radius[repetition], effect[repetition] + reference_radius[repetition])
        if "reference_point" in data:
            reference = (data["reference_point"][repetition]-data["reference_width"][repetition]/2,
                         data["reference_point"][repetition]+data["reference_width"][repetition]/2)
        result = protected_projection(mean[repetition], region, residual_sets, moments,
                                      prediction_radius, reference, fallback_rule)
        intervals[repetition], points[repetition], codes[repetition], envelopes[repetition] = result
    if "parameter_bounds" in data:
        intervals = np.clip(intervals,*data["parameter_bounds"])
        points = np.minimum(intervals[:,1],np.maximum(intervals[:,0],points))
    return dict(interval=intervals, point_midpoint=intervals.mean(axis=1), point_source=points,
                fallback=codes, length_envelope=envelopes,
                alpha_target=target_alpha, alpha_source_arm=alpha_arm,
                source_bias_budget=source_budget, target_bias_budget=target_budget)

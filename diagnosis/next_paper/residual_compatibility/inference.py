"""Research-only residual confidence sets and exact ellipsoid projection.

No fitted RoCE inference claim. Source mean noise assumptions and calibration
modes are explicit; source arm sets may differ and target covariance is joint.
"""
import math
import numpy as np
from scipy.stats import chi2, norm, ncx2


def intersect_sets(first, second, tolerance=1e-12):
    result = []
    i = j = 0
    while i < len(first) and j < len(second):
        lower = max(first[i][0], second[j][0])
        upper = min(first[i][1], second[j][1])
        if lower <= upper + tolerance:
            result.append((lower - tolerance, max(lower, upper) + tolerance))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return result


def fourier_interval_set(transform, frequencies, allowance, valid_fraction, domain):
    """Invert every complex disk as analytic arcs; never use a search grid."""
    if not 0 < valid_fraction <= 1 or np.any(np.asarray(frequencies) <= 0):
        raise ValueError("Positive frequencies and a validity fraction in (0,1] required")
    accepted = [tuple(domain)]
    for value, frequency, slack in zip(transform, frequencies, allowance):
        magnitude = abs(value)
        radius = 1 - valid_fraction + slack
        if magnitude < 1e-15:
            if valid_fraction > radius + 1e-12:
                return []
            continue
        threshold = (magnitude**2 + valid_fraction**2 - radius**2) / (2 * magnitude * valid_fraction)
        if threshold <= -1:
            continue
        if threshold > 1 + 1e-12:
            return []
        half_angle = math.acos(min(1, threshold))
        phase = np.angle(value)
        first_period = math.ceil((frequency * domain[0] - phase - half_angle) / (2 * math.pi))
        last_period = math.floor((frequency * domain[1] - phase + half_angle) / (2 * math.pi))
        arcs = [((phase + 2 * math.pi * period - half_angle) / frequency,
                 (phase + 2 * math.pi * period + half_angle) / frequency)
                for period in range(first_period, last_period + 1)]
        accepted = intersect_sets(accepted, arcs)
        if not accepted:
            return []
    return accepted


def fourier_summary(estimates, variances, alpha, frequency_scale=2.0,
                    calibration="bernstein", variance_upper=None, score_ranges=None,
                    sample_sizes=None):
    """Batched characteristic functions. variances are fixed before evaluation.

    With variance_upper, variances must be certified LOWER endpoints. The
    resulting attenuation error is added explicitly. Bounded-score mode adds
    a Lindeberg characteristic-function remainder, and uses a distribution-free
    variance bound for complex exponentials. Missing remainder is not inferred.
    """
    estimates = np.asarray(estimates)
    if estimates.ndim != 2 or not np.all(np.isfinite(estimates)) or not 0 < alpha < 1:
        raise ValueError("Finite repeat-by-source estimates and alpha in (0,1) required")
    if calibration not in ["hoeffding", "bernstein"] or frequency_scale <= 0:
        raise ValueError("Unknown calibration or nonpositive frequency scale")
    variances = np.broadcast_to(variances, estimates.shape)
    if not np.all(np.isfinite(variances)):
        raise ValueError("Finite variances required")
    if variance_upper is not None:
        variance_upper = np.broadcast_to(variance_upper, estimates.shape)
        if np.any(variance_upper < variances) or not np.all(np.isfinite(variance_upper)):
            raise ValueError("Variance upper endpoints must dominate lower endpoints")
    if score_ranges is not None and (sample_sizes is None or np.any(np.asarray(sample_sizes) < 2)
                                     or np.any(np.asarray(score_ranges) < 0)):
        raise ValueError("Bounded-score calibration requires ranges and at least two patients/site")
    count = estimates.shape[1]
    maximum = np.max(variances if variance_upper is None else variance_upper, axis=1)
    if np.any(maximum <= 0) or np.any(variances < 0):
        raise ValueError("Positive variance scale and nonnegative variances required")
    levels = np.arange(1, int(math.ceil(math.log(math.e * count))) + 1)
    frequencies = np.sqrt(levels)[None, :] / (frequency_scale * np.sqrt(maximum)[:, None])
    transforms = np.empty_like(frequencies, dtype=complex)
    slacks = np.empty_like(frequencies)
    logarithm = math.log(4 * len(levels) / alpha)
    for column in range(len(levels)):
        frequency = frequencies[:, column, None]
        amplitude = np.exp(frequency**2 * variances / 2)
        transforms[:, column] = np.mean(amplitude * np.exp(1j * frequency * estimates), axis=1)
        hoeffding = 2 * np.sqrt(logarithm * np.sum(amplitude**2, axis=1)) / count
        if score_ranges is not None:
            noise_variance_bound = variances if variance_upper is None else variance_upper
            complex_variance = amplitude**2 * np.minimum(1, frequency**2 * noise_variance_bound)
            bound = 2 * np.max(amplitude, axis=1)
        elif variance_upper is not None:
            complex_variance = amplitude**2 * (-np.expm1(-frequency**2 * variance_upper))
            bound = np.max(amplitude, axis=1) + 1
        else:
            complex_variance = np.expm1(frequency**2 * variances)
            bound = np.max(amplitude, axis=1) + 1
        term = bound * logarithm / 3
        bernstein = math.sqrt(2) * (term + np.sqrt(term**2 +
            2 * np.sum(complex_variance, axis=1) * logarithm)) / count
        slack = hoeffding if calibration == "hoeffding" else np.minimum(hoeffding, bernstein)
        if variance_upper is not None:
            slack += np.mean(-np.expm1(-frequency**2 * (variance_upper - variances) / 2), axis=1)
        if score_ranges is not None:
            sizes = np.broadcast_to(sample_sizes, estimates.shape)
            patient_variance = (variances if variance_upper is None else variance_upper) * sizes
            third_moment_bound = score_ranges * patient_variance + \
                2 * math.sqrt(2 / math.pi) * patient_variance**1.5
            approximation = amplitude * abs(frequency)**3 * third_moment_bound / (6 * sizes**2)
            slack += np.mean(approximation, axis=1)
        slacks[:, column] = slack
    return transforms, frequencies, slacks


def count_interval_set(estimates, variances, valid_minimum, alpha, domain):
    radius = norm.isf(alpha / (2 * len(estimates))) * np.sqrt(variances)
    events = []
    for lower, upper in zip(estimates - radius, estimates + radius):
        lower, upper = max(lower, domain[0]), min(upper, domain[1])
        if lower <= upper:
            events.extend([(lower, 1), (upper, -1)])
    events.sort(key=lambda item: (item[0], -item[1]))
    active = 0
    start = None
    result = []
    for position, change in events:
        previous = active
        active += change
        if previous < valid_minimum <= active:
            start = position
        if active < valid_minimum <= previous:
            result.append((start, position))
    return result


def noncentrality_upper(statistics, degrees, alpha):
    statistics = np.maximum(statistics, 0)
    lower = np.zeros_like(statistics)
    upper = np.maximum(1, statistics - degrees + 10 * np.sqrt(statistics + degrees) + 20)
    active = ncx2.cdf(statistics, degrees, 0) > alpha
    for _ in range(40):
        enlarge = active & (ncx2.cdf(statistics, degrees, upper) > alpha)
        if not np.any(enlarge):
            break
        upper[enlarge] *= 2
    else:
        raise ArithmeticError("Could not bracket noncentrality")
    for _ in range(45):
        middle = (lower + upper) / 2
        above = ncx2.cdf(statistics, degrees, middle) > alpha
        lower = np.where(active & above, middle, lower)
        upper = np.where(active & ~above, middle, upper)
    return np.where(active, upper, 0)


def dispersion_intervals(estimates, variances, valid_minimum, alpha,
                         sample_variances=None, sample_sizes=None, score_ranges=None,
                         variance_status="known_gaussian"):
    """Known Gaussian variances or certified bounds with unbiased correction.

    The exact noncentral chi-square branch is reserved for deterministic known
    Gaussian variances. A random variance subtraction is a different statistic.
    gaussian_plugin is diagnostic only and carries no coverage assertion.
    """
    if variance_status not in ["known_gaussian", "certified", "gaussian_plugin"]:
        raise ValueError("Unknown variance status")
    if variance_status == "certified" and sample_variances is None:
        raise ValueError("Certified variance bounds require unbiased evaluation correction")
    count = estimates.shape[1]
    variances = np.broadcast_to(variances, estimates.shape)
    center = np.mean(estimates, axis=1)
    variance_mean = np.sum(variances, axis=1) / count**2
    if score_ranges is None:
        noise_radius = norm.isf(alpha / 4) * np.sqrt(variance_mean)
    else:
        logarithm = math.log(4 / alpha)
        patient_ranges = np.broadcast_to(score_ranges, estimates.shape)
        sizes = np.broadcast_to(sample_sizes, estimates.shape)
        bound = np.max(patient_ranges / (count * sizes), axis=1)
        term = bound * logarithm / 3
        noise_radius = term + np.sqrt(term**2 + 2 * variance_mean * logarithm)
    if valid_minimum == count:
        return np.column_stack((center - noise_radius, center + noise_radius))
    observed_variance = np.mean((estimates - center[:, None])**2, axis=1)
    homogeneous = np.max(np.ptp(variances, axis=1)) < 1e-14
    if homogeneous and sample_variances is None and score_ranges is None and variance_status == "known_gaussian":
        variance = variances[:, 0]
        statistic = count * observed_variance / variance
        energy_upper = variance / count * noncentrality_upper(statistic, count - 1, alpha / 2)
    else:
        correction = variances if sample_variances is None else sample_variances
        observed = observed_variance - (count - 1) / count**2 * np.sum(correction, axis=1)
        kappa = 1 if sample_variances is None else sample_sizes / (sample_sizes - 1)
        linear = 4 * np.max(variances, axis=1) / count
        quadratic = 2 / count**4 * (np.sum(variances, axis=1)**2 +
            np.sum((kappa * (count - 1)**2 - 1) * variances**2, axis=1))
        multiplier = (1 - alpha / 2) / (alpha / 2)
        energy_upper = np.maximum(0, observed + multiplier * linear / 2 +
            np.sqrt(np.maximum(0, multiplier * linear * observed +
                               multiplier**2 * linear**2 / 4 + multiplier * quadratic)))
    radius = noise_radius + np.sqrt((count - valid_minimum) / valid_minimum * energy_upper)
    return np.column_stack((center - radius, center + radius))


class TargetRegion:
    """Ellipsoid for (prediction contrast, treated residual, control residual).

    Projection over a rectangle needs only the unconstrained optimum and four
    edges. Every edge is solved analytically by conditioning a Gaussian ellipse.
    The same geometry applies to moment ellipsoids. Singular known prediction
    contrasts are allowed; both residual variances must be positive definite.
    """
    def __init__(self, covariance, alpha=.025, calibration="gaussian"):
        if calibration not in ["gaussian", "moment"] or not 0 < alpha < 1:
            raise ValueError("Unknown target calibration or invalid alpha")
        self.covariance = np.asarray(covariance, dtype=float)
        eigenvalues = np.linalg.eigvalsh(self.covariance)
        if np.min(eigenvalues) < -1e-12 or np.min(np.linalg.eigvalsh(self.covariance[1:, 1:])) <= 0:
            raise ValueError("PSD covariance and positive residual block required")
        rank = np.sum(eigenvalues > np.max(eigenvalues) * 1e-10)
        self.critical_squared = chi2.isf(alpha, rank) if calibration == "gaussian" else rank / alpha
        self.contrast = np.array([1., 1., -1.])
        self.effect_variance = float(self.contrast @ self.covariance @ self.contrast)
        self.offset = math.sqrt(self.critical_squared / self.effect_variance) * self.covariance @ self.contrast
        self.edges = []
        for fixed in [1, 2]:
            other = 3 - fixed
            remaining = [0, other]
            variance = self.covariance[fixed, fixed]
            slope = self.covariance[remaining, fixed] / variance
            conditional = self.covariance[np.ix_(remaining, remaining)] - np.outer(slope, slope) * variance
            contrast = self.contrast[remaining]
            contrast_variance = max(0, float(contrast @ conditional @ contrast))
            self.edges.append((fixed, other, remaining, variance, slope, conditional, contrast, contrast_variance))

    def domains(self, center):
        radii = np.sqrt(self.critical_squared * np.diag(self.covariance)[1:])
        return [(center[index + 1] - radii[index], center[index + 1] + radii[index]) for index in range(2)]

    def project(self, center, residual_sets):
        center = np.asarray(center)
        values = []
        if not residual_sets[0] or not residual_sets[1]:
            return None
        for treated in residual_sets[0]:
            for control in residual_sets[1]:
                rectangle = [None, treated, control]
                for sign in [-1, 1]:
                    optimum = center + sign * self.offset
                    if all(rectangle[index][0] <= optimum[index] <= rectangle[index][1] for index in [1, 2]):
                        values.append(float(self.contrast @ optimum))
                for fixed, other, remaining, variance, slope, conditional, contrast, contrast_variance in self.edges:
                    for boundary in rectangle[fixed]:
                        displacement = boundary - center[fixed]
                        budget = self.critical_squared - displacement**2 / variance
                        if budget < -1e-10:
                            continue
                        budget = max(0, budget)
                        conditional_mean = center[remaining] + slope * displacement
                        other_radius = math.sqrt(max(0, budget * conditional[1, 1]))
                        lower = max(rectangle[other][0], conditional_mean[1] - other_radius)
                        upper = min(rectangle[other][1], conditional_mean[1] + other_radius)
                        if lower > upper + 1e-12:
                            continue
                        upper = max(lower, upper)
                        for sign in [-1, 1]:
                            other_optimum = conditional_mean[1]
                            if contrast_variance > 1e-20:
                                other_optimum += sign * math.sqrt(budget / contrast_variance) * (conditional @ contrast)[1]
                            other_value = min(upper, max(lower, other_optimum))
                            remaining_budget = max(0, budget - (other_value - conditional_mean[1])**2 / conditional[1, 1])
                            prediction_mean = conditional_mean[0] + conditional[0, 1] / conditional[1, 1] * (other_value - conditional_mean[1])
                            prediction_variance = max(0, conditional[0, 0] - conditional[0, 1]**2 / conditional[1, 1])
                            prediction_value = prediction_mean + sign * math.sqrt(prediction_variance * remaining_budget)
                            values.append(prediction_value + self.contrast[fixed] * boundary + self.contrast[other] * other_value)
        if not values:
            return None
        return min(values) - 1e-11, max(values) + 1e-11

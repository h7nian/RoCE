"""Finite-stratum nuisance certificates for unknown valid source subsets.

Training determines weights and predictions. Held-out source X,A data
certify linear calibration moments at all outcome-error-box corners. The target
probability region is simultaneous, so its fitted coefficients may be random.
No assumption is made that every source is outcome-transportable.
"""
import itertools
import numpy as np
from concentrated_bounds import (pooled_mean_interval, exponential_dispersion_upper,
                                 simplex_linear_max, grouped_probability_max)
from fitted_strata import binomial_bounds


def calibration_remainder_budgets(data, valid_minimum, failure=.005):
    if not 0 < failure < 1:
        raise ValueError("Positive failure probability below one required")
    weights = data["source_weights"]
    repeats, count, cells, arms = weights.shape
    if arms != 2 or cells > 8:
        raise ValueError("This exact corner diagnostic supports two arms and at most eight strata")
    valid = np.broadcast_to(valid_minimum, (2,))
    if np.any(valid < 1) or np.any(valid > count) or np.any(valid != np.floor(valid)):
        raise ValueError("Invalid guaranteed source count")
    training_counts = data["target_training_counts"]
    lower, upper = binomial_bounds(training_counts[..., 1], training_counts.sum(axis=3),
                                   failure * .2 / (2 * cells))
    lower -= data["common_prediction"]
    upper -= data["common_prediction"]
    joint_counts = data["target_full_joint_counts"]
    total = joint_counts.sum(axis=(1, 2))
    joint_lower, joint_upper = binomial_bounds(joint_counts, total[:, None, None],
                                               failure * .1 / (2 * cells))
    cell_lower, cell_upper = binomial_bounds(joint_counts.sum(axis=2), total[:, None],
                                             failure * .1 / cells)
    effective_lower = np.maximum(cell_lower, joint_lower.sum(axis=2))
    effective_upper = np.minimum(cell_upper, joint_upper.sum(axis=2))
    corners = list(itertools.product([0, 1], repeat=cells))
    corner_failure = failure * .6 / (2 * len(corners))
    source_mean_budget = np.zeros((repeats, 2))
    source_squared_budget = np.zeros((repeats, 2))
    target_budget = np.zeros((repeats, 3))
    frequencies = data.get("source_design_joint_frequencies", data["pilot_joint_frequencies"])
    design_size = data.get("source_design_sample_size", data["pilot_size"])
    for arm in [0, 1]:
        for corner in corners:
            error = np.where(np.array(corner)[None, :], upper[:, :, arm], lower[:, :, arm])
            score = weights[:, :, :, arm] * error[:, None, :]
            design_mean = np.sum(frequencies[:, :, :, arm] * score, axis=2)
            design_variance = (np.sum(frequencies[:, :, :, arm] * score**2, axis=2) - design_mean**2) / (design_size - 1)
            score_range = np.maximum(0, score.max(axis=2)) - np.minimum(0, score.min(axis=2))
            interval = pooled_mean_interval(design_mean, design_variance, design_size,
                                            np.max(abs(score), axis=2), corner_failure / 2)
            energy = exponential_dispersion_upper(design_mean, design_variance, design_size,
                                                  score_range, corner_failure / 2)
            target_upper = simplex_linear_max(error, effective_lower, effective_upper)
            target_lower = -simplex_linear_max(-error, effective_lower, effective_upper)
            average = np.maximum(abs(interval[:, 0] - target_upper), abs(interval[:, 1] - target_lower))
            # This protects every subset of g sources, including the unknown
            # truly valid subset, rather than assuming its average is observable.
            source_mean_budget[:, arm] = np.maximum(source_mean_budget[:, arm],
                average + np.sqrt((count - valid[arm]) / valid[arm] * energy))
            source_squared_budget[:, arm] = np.maximum(source_squared_budget[:, arm], energy + average**2)
            coefficients = np.repeat(-error[:, :, None], 2, axis=2)
            coefficients[:, :, arm] += error * data["target_inverse"][:, :, arm]
            for sign in [-1, 1]:
                bound = grouped_probability_max(sign * coefficients, joint_lower, joint_upper,
                                                 cell_lower, cell_upper)
                target_budget[:, arm + 1] = np.maximum(target_budget[:, arm + 1], bound)
    return dict(source_mean_budget=source_mean_budget, source_squared_budget=source_squared_budget,
                target_budget=target_budget, failure=failure, corner_count=len(corners))


def target_remainder_budget(training_counts, full_joint_counts, failure):
    """Joint-probability certificate for one saturated target fit, both arms."""
    arm_counts = training_counts.sum(axis=3)
    cell_counts = arm_counts.sum(axis=2)
    common = (training_counts[..., 1] + .5) / (arm_counts + 1)
    assignment = (arm_counts[:, :, 0] + .5) / (cell_counts + 1)
    inverse = 1 / np.stack((assignment, 1 - assignment), axis=2)
    repeats, cells = cell_counts.shape
    lower, upper = binomial_bounds(training_counts[..., 1], arm_counts, failure / (4 * cells))
    lower -= common
    upper -= common
    total = full_joint_counts.sum(axis=(1, 2))
    joint_lower, joint_upper = binomial_bounds(full_joint_counts, total[:, None, None], failure / (8 * cells))
    cell_lower, cell_upper = binomial_bounds(full_joint_counts.sum(axis=2), total[:, None], failure / (4 * cells))
    budget = np.zeros((repeats, 2))
    for arm in [0, 1]:
        for corner in itertools.product([0, 1], repeat=cells):
            error = np.where(np.array(corner)[None, :], upper[:, :, arm], lower[:, :, arm])
            coefficients = np.repeat(-error[:, :, None], 2, axis=2)
            coefficients[:, :, arm] += error * inverse[:, :, arm]
            for sign in [-1, 1]:
                budget[:, arm] = np.maximum(budget[:, arm], grouped_probability_max(
                    sign * coefficients, joint_lower, joint_upper, cell_lower, cell_upper))
    return budget

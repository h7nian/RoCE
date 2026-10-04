"""Finite-binomial certificates for a fixed but unknown valid source subset.

This diagnostic keeps the existing half-count weights. It requires a declared
lower bound on true source joint (cell,arm) probabilities. No validity set is
selected from fitted errors. Source training is averaged before bounding its
nuisance contribution; the target training component remains shared.
"""
import itertools
import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom
from fitted_strata import binomial_bounds
from concentrated_bounds import simplex_linear_max, grouped_probability_max


def binomial_ratio_mgf(training_size, cell_count, probability_minimum,
                       grid_size=2048, tilts=None):
    """Uniform MGF envelope via stochastic bracketing, not pointwise p gridding.

    For p in [l,u], B(l) <=st B(p) <=st B(u). Therefore (n+L)p/(B(p)+.5)
    is bracketed by (n+L)l/(B(u)+.5) and (n+L)u/(B(l)+.5). Endpoints give
    conservative finite-grid MGF envelopes for every probability in the bin.
    """
    if (training_size < 1 or training_size != int(training_size)
            or cell_count < 1 or cell_count != int(cell_count)
            or not 0 < probability_minimum < 1 or grid_size < 2):
        raise ValueError("Positive integer sizes and probability minimum in (0,1) required")
    if tilts is None:
        tilts = np.geomspace(.001, 2., 100)
    tilts = np.asarray(tilts, dtype=float)
    if np.any(tilts <= 0) or not np.all(np.isfinite(tilts)):
        raise ValueError("Finite positive deterministic tilts required")
    # Logarithmic probability cells give comparable relative precision near p_min.
    edges = np.geomspace(probability_minimum, 1., grid_size + 1)
    lower, upper = edges[:-1, None], edges[1:, None]
    successes = np.arange(training_size + 1)[None, :]
    lower_log_probability = binom.logpmf(successes, training_size, lower)
    upper_log_probability = binom.logpmf(successes, training_size, upper)
    ratio_upper = (training_size + cell_count) * upper / (successes + .5)
    ratio_lower = (training_size + cell_count) * lower / (successes + .5)
    log_upper = []
    log_lower = []
    for tilt in tilts:
        log_upper.append(np.max(logsumexp(lower_log_probability + tilt * (ratio_upper - 1), axis=1)))
        log_lower.append(np.max(logsumexp(upper_log_probability - tilt * (ratio_lower - 1), axis=1)))
    # Floating-point cushion is applied outwards. This is analytic stochastic
    # bracketing evaluated numerically, not interval-arithmetic certification.
    return dict(tilts=tilts, log_upper=np.array(log_upper) + 1e-10,
                log_lower=np.array(log_lower) + 1e-10,
                training_size=training_size, cell_count=cell_count,
                probability_minimum=probability_minimum, grid_size=grid_size)


def ratio_average_interval(envelope, valid_count, failure):
    """Two-sided bound for an unobserved fixed g-subset's average ratio."""
    if valid_count < 1 or valid_count != int(valid_count) or not 0 < failure < 1:
        raise ValueError("Positive integer subset size and failure in (0,1) required")
    logarithm = np.log(2 / failure) / valid_count
    upper = np.min((envelope['log_upper'] + logarithm) / envelope['tilts'])
    lower = -np.min((envelope['log_lower'] + logarithm) / envelope['tilts'])
    return float(1 + lower), float(1 + upper)


def training_average_budgets(data, valid_minimum, envelope, failure=.005):
    if data.get('nuisance_structure') != 'saturated_design_weights':
        raise ValueError('This bound requires saturated weights fitted without outcomes')
    weights = data['source_weights']
    repeats, count, cells, arms = weights.shape
    if arms != 2 or cells > 8 or envelope['cell_count'] != cells:
        raise ValueError("Envelope must match the two-arm finite-stratum model")
    if envelope['training_size'] != data['source_training_size']:
        raise ValueError("Envelope and source training sizes differ")
    valid = np.broadcast_to(valid_minimum, (2,))
    if np.any(valid < 1) or np.any(valid > count) or np.any(valid != np.floor(valid)):
        raise ValueError("Invalid fixed valid subset size")
    training = data['target_training_counts']
    arm_counts = training.sum(axis=3)
    cell_counts = arm_counts.sum(axis=2)
    fitted_mass = (cell_counts + .5) / (data['target_training_size'] + .5 * cells)
    lower, upper = binomial_bounds(training[..., 1], arm_counts, failure * .2 / (2 * cells))
    lower -= data['common_prediction']
    upper -= data['common_prediction']
    joint_counts = data['target_full_joint_counts']
    total = joint_counts.sum(axis=(1, 2))
    joint_lower, joint_upper = binomial_bounds(joint_counts, total[:, None, None], failure * .1 / (2 * cells))
    cell_lower, cell_upper = binomial_bounds(joint_counts.sum(axis=2), total[:, None], failure * .1 / cells)
    effective_lower = np.maximum(cell_lower, joint_lower.sum(axis=2))
    effective_upper = np.minimum(cell_upper, joint_upper.sum(axis=2))
    source_budget = np.zeros((repeats, 2))
    target_budget = np.zeros((repeats, 3))
    ratio_bounds = np.array([ratio_average_interval(envelope, int(g), failure * .6 / (2 * cells)) for g in valid])
    shared_budget = np.zeros_like(source_budget)
    for arm in [0, 1]:
        for corner in itertools.product([0, 1], repeat=cells):
            error = np.where(np.array(corner)[None, :], upper[:, :, arm], lower[:, :, arm])
            target_upper = simplex_linear_max(error, effective_lower, effective_upper)
            target_lower = -simplex_linear_max(-error, effective_lower, effective_upper)
            source_upper = np.sum(fitted_mass * error * np.where(error >= 0, ratio_bounds[arm, 1], ratio_bounds[arm, 0]), axis=1)
            source_lower = np.sum(fitted_mass * error * np.where(error >= 0, ratio_bounds[arm, 0], ratio_bounds[arm, 1]), axis=1)
            source_budget[:, arm] = np.maximum(source_budget[:, arm], np.maximum(source_upper - target_lower, target_upper - source_lower))
            shared = np.sum(fitted_mass * error, axis=1)
            shared_budget[:, arm] = np.maximum(shared_budget[:, arm], np.maximum(shared - target_lower, target_upper - shared))
            coefficients = np.repeat(-error[:, :, None], 2, axis=2)
            coefficients[:, :, arm] += error * data['target_inverse'][:, :, arm]
            for sign in [-1, 1]:
                target_budget[:, arm + 1] = np.maximum(target_budget[:, arm + 1], grouped_probability_max(
                    sign * coefficients, joint_lower, joint_upper, cell_lower, cell_upper))
    return dict(source_mean_budget=source_budget, target_budget=target_budget,
                ratio_bounds=ratio_bounds, shared_target_budget=shared_budget)

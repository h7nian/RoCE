"""Scalar nuisance drift bounds conditional on outcome-independent designs.

Valid only when the fitted inverse weights depend on X,A, not target training
outcomes. Current saturated half-count fits satisfy this restriction. The true
valid source subset is fixed before sampling; it is never chosen by the code.
"""
import numpy as np
from fitted_strata import binomial_bounds
from concentrated_bounds import pooled_mean_interval, exponential_dispersion_upper
from training_average import ratio_average_interval


def outcome_linear_radius(coefficient_bound, arm_counts, failure):
    """Hoeffding for a single unknown but design-fixed coefficient vector.

    Half-count fitted cell means have conditional smoothing bias bounded by
    1/[2(N+1)]. The noise has proxy sum N*c^2/[4(N+1)^2]. No union over
    outcome-box corners is required; its direction is fixed given the design.
    """
    bound, counts = np.broadcast_arrays(coefficient_bound, arm_counts)
    if (not 0 < failure < 1 or np.any(bound < 0) or np.any(counts < 0)
            or np.any(counts != np.floor(counts)) or not np.all(np.isfinite(bound)) or not np.all(np.isfinite(counts))):
        raise ValueError("Finite nonnegative coefficients/counts and failure in (0,1) required")
    smoothing = .5 * np.sum(bound / (counts + 1), axis=-1)
    noise = np.sqrt(.5 * np.log(2 / failure) * np.sum(counts * bound**2 / (counts + 1)**2, axis=-1))
    return smoothing + noise


def conditional_design_coefficients(data, valid_minimum, target_failure, source_failure,
                                    source_calibration='empirical_design', envelope=None):
    if source_calibration not in ['empirical_design', 'training_mgf'] or not 0 < target_failure < 1 or not 0 < source_failure < 1:
        raise ValueError("Unknown source calibration or invalid allowance")
    if data.get('nuisance_structure') != 'saturated_design_weights':
        raise ValueError('This bound requires saturated weights fitted without outcomes')
    weights = data['source_weights']
    repeats, count, cells, arms = weights.shape
    valid = np.broadcast_to(valid_minimum, (2,))
    if arms != 2 or np.any(valid < 1) or np.any(valid > count) or np.any(valid != np.floor(valid)):
        raise ValueError("Two arms and valid subset sizes in [1,K] required")
    if source_calibration == 'training_mgf':
        if envelope is None or envelope['cell_count'] != cells or envelope['training_size'] != data['source_training_size']:
            raise ValueError("Matching training MGF envelope required")
    counts = data['target_training_counts'].sum(axis=3)
    joint_counts = data['target_full_joint_counts']
    total = joint_counts.sum(axis=(1, 2))
    joint_lower, joint_upper = binomial_bounds(joint_counts, total[:, None, None], target_failure / (4 * cells))
    cell_lower, cell_upper = binomial_bounds(joint_counts.sum(axis=2), total[:, None], target_failure / (2 * cells))
    cell_lower = np.maximum(cell_lower, joint_lower.sum(axis=2))
    cell_upper = np.minimum(cell_upper, joint_upper.sum(axis=2))
    fitted_mass = (counts.sum(axis=2) + .5) / (data['target_training_size'] + .5 * cells)
    source_coefficients = np.empty((repeats, cells, 2))
    target_coefficients = np.empty_like(source_coefficients)
    cell_failure = source_failure / (2 * cells)
    design_size = data.get('source_design_sample_size', data['pilot_size'])
    design_frequencies = data.get('source_design_joint_frequencies', data['pilot_joint_frequencies'])
    for arm in [0, 1]:
        for cell in range(cells):
            if source_calibration == 'empirical_design':
                score = weights[:, :, cell, arm]
                frequency = design_frequencies[:, :, cell, arm]
                mean = score * frequency
                variance = score**2 * frequency * (1 - frequency) / (design_size - 1)
                interval = pooled_mean_interval(mean, variance, design_size, score, cell_failure / 2)
                energy = exponential_dispersion_upper(mean, variance, design_size, score, cell_failure / 2)
                radius = np.sqrt((count - valid[arm]) / valid[arm] * energy)
                source_lower, source_upper = interval[:, 0] - radius, interval[:, 1] + radius
            else:
                ratio_lower, ratio_upper = ratio_average_interval(envelope, int(valid[arm]), cell_failure)
                source_lower = fitted_mass[:, cell] * ratio_lower
                source_upper = fitted_mass[:, cell] * ratio_upper
            source_coefficients[:, cell, arm] = np.maximum(abs(source_lower - cell_upper[:, cell]),
                                                          abs(source_upper - cell_lower[:, cell]))
        inverse = data['target_inverse'][:, :, arm]
        other = 1 - arm
        # Exact two-variable linear extrema, respecting the cell marginal band.
        first_high = np.minimum(joint_upper[:, :, arm], cell_upper - joint_lower[:, :, other])
        second_low = np.maximum(joint_lower[:, :, other], cell_lower - first_high)
        first_low = np.maximum(joint_lower[:, :, arm], cell_lower - joint_upper[:, :, other])
        second_high = np.minimum(joint_upper[:, :, other], cell_upper - first_low)
        coefficient_upper = (inverse - 1) * first_high - second_low
        coefficient_lower = (inverse - 1) * first_low - second_high
        target_coefficients[:, :, arm] = np.maximum(abs(coefficient_lower), abs(coefficient_upper))
    return dict(source_coefficient_bound=source_coefficients, target_coefficient_bound=target_coefficients,
                target_design_failure=target_failure, source_design_failure=source_failure,
                valid_minimum=valid.copy(), source_calibration=source_calibration)


def conditional_remainder_budgets(data, valid_minimum, failure=.005,
                                  source_calibration='empirical_design', envelope=None):
    if not 0 < failure < 1:
        raise ValueError("Confidence failure must be in (0,1)")
    coefficients = conditional_design_coefficients(data, valid_minimum, failure * .2,
        failure * .4, source_calibration, envelope)
    source_coefficients = coefficients['source_coefficient_bound']
    target_coefficients = coefficients['target_coefficient_bound']
    counts = data['target_training_counts'].sum(axis=3)
    repeats = len(counts)
    # Four scalar outcome events: source subset and target, in each arm.
    scalar_failure = failure * .4 / 4
    source_budget = np.column_stack([outcome_linear_radius(source_coefficients[:, :, arm], counts[:, :, arm], scalar_failure)
                                     for arm in [0, 1]])
    target_budget = np.column_stack([np.zeros(repeats)] + [outcome_linear_radius(
        target_coefficients[:, :, arm], counts[:, :, arm], scalar_failure) for arm in [0, 1]])
    return dict(source_mean_budget=source_budget, target_budget=target_budget,
                source_coefficient_bound=source_coefficients, target_coefficient_bound=target_coefficients,
                failure=failure, source_calibration=source_calibration)


def conditional_target_budget(training_counts, full_joint_counts, failure):
    """Same scalar improvement for one target-only fold, without any sources."""
    counts = training_counts.sum(axis=3)
    repeats, cells, arms = counts.shape
    total = full_joint_counts.sum(axis=(1, 2))
    lower, upper = binomial_bounds(full_joint_counts, total[:, None, None], failure * .25 / (2 * cells))
    mass_lower, mass_upper = binomial_bounds(full_joint_counts.sum(axis=2), total[:, None], failure * .25 / cells)
    mass_lower = np.maximum(mass_lower, lower.sum(axis=2))
    mass_upper = np.minimum(mass_upper, upper.sum(axis=2))
    propensity = (counts[:, :, 0] + .5) / (counts.sum(axis=2) + 1)
    inverse = 1 / np.stack((propensity, 1 - propensity), axis=2)
    budgets = np.zeros((repeats, 2))
    for arm in [0, 1]:
        other = 1 - arm
        first_high = np.minimum(upper[:, :, arm], mass_upper - lower[:, :, other])
        second_low = np.maximum(lower[:, :, other], mass_lower - first_high)
        first_low = np.maximum(lower[:, :, arm], mass_lower - upper[:, :, other])
        second_high = np.minimum(upper[:, :, other], mass_upper - first_low)
        bound = np.maximum(abs((inverse[:, :, arm] - 1) * first_high - second_low),
                           abs((inverse[:, :, arm] - 1) * first_low - second_high))
        budgets[:, arm] = outcome_linear_radius(bound, counts[:, :, arm], failure * .25)
    return budgets


def full_target_conditional_interval(data):
    """Fair 1000-patient two-fold target reference using the same improvement."""
    if data.get('nuisance_structure') != 'saturated_design_weights':
        raise ValueError('Conditional target certificate requires design-only weights')
    halves = [data['target_training_counts'], data['target_evaluation_counts']]
    cells = halves[0].shape[1]
    states = np.array([(cell, arm, outcome) for cell in range(cells) for arm in [0, 1] for outcome in [0, 1]])
    index, arm_index, outcome = states.T
    radii = []
    means = []
    for train, evaluate in [halves, halves[::-1]]:
        arm_counts = train.sum(axis=3)
        common = (train[..., 1] + .5) / (arm_counts + 1)
        propensity = (arm_counts[:, :, 0] + .5) / (arm_counts.sum(axis=2) + 1)
        inverse = 1 / np.stack((propensity, 1 - propensity), axis=2)
        score = common[:, index, 0] - common[:, index, 1]
        for arm, sign in [(0, 1), (1, -1)]:
            score += sign * (arm_index == arm)[None, :] * inverse[:, index, arm] * (outcome[None, :] - common[:, index, arm])
        counts = evaluate.reshape(len(evaluate), -1)
        size = counts.sum(axis=1)
        mean = np.sum(counts * score, axis=1) / size
        variance = (np.sum(counts * score**2, axis=1) - size * mean**2) / (size * (size - 1))
        logarithm = np.log(4 / .0225)
        radius = np.sqrt(2 * np.maximum(variance, 0) * logarithm) + 7 * np.ptp(score, axis=1) * logarithm / (3 * (size - 1))
        radius += conditional_target_budget(train, data['target_full_joint_counts'], .0025).sum(axis=1)
        radii.append(radius)
        means.append(mean)
    point = np.mean(means, axis=0)
    if not np.allclose(point, data['reference_point'], atol=1e-12, rtol=1e-12):
        raise ArithmeticError('Reconstructed target folds differ from the reference point')
    radius = np.mean(radii, axis=0)
    return np.clip(np.column_stack((point - radius, point + radius)), -1, 1)


def reuse_source_design_for_evaluation(data):
    """Merge the two held-out source halves for direct dispersion only.

    Nuisance certificates and dispersion events are combined by a union bound,
    so they may reuse observations. No evaluation-dependent Fourier frequency
    or variance certificate is used. The initial fits stay on their 500 patients.
    """
    pilot = data['pilot_size']
    evaluation = data['source_evaluation_size']
    if pilot < 2 or evaluation < 2:
        raise ValueError('Two held-out source blocks with at least two patients required')
    total = pilot + evaluation
    mean = (pilot * data['pilot_source_mean'] + evaluation * data['source']) / total
    # pilot_variance is S_pilot^2 / old evaluation size in the legacy packets.
    within = (pilot - 1) * evaluation * data['pilot_variance'] + (evaluation - 1) * evaluation * data['sample_variance']
    between = pilot * (data['pilot_source_mean'] - mean)**2 + evaluation * (data['source'] - mean)**2
    variance = (within + between) / (total * (total - 1))
    frequencies = (pilot * data['pilot_joint_frequencies'] + evaluation * data['evaluation_joint_frequencies']) / total
    extra = {}
    if 'source_sample_cross_covariance' in data:
        cross_within = ((pilot-1)*evaluation*data['pilot_sample_cross_covariance']
                        + (evaluation-1)*evaluation*data['source_sample_cross_covariance'])
        cross_between = (pilot*np.prod(data['pilot_source_mean']-mean,axis=2)
                         + evaluation*np.prod(data['source']-mean,axis=2))
        extra['source_sample_cross_covariance'] = (cross_within+cross_between)/(total*(total-1))
    return dict(data, **extra, source=mean, sample_variance=variance,
        variance=data['variance'] * evaluation / total,
        sample_sizes=np.full(data['source'].shape[1], total), source_evaluation_size=total,
        source_design_sample_size=total, source_design_joint_frequencies=frequencies,
        pilot_size=0, source_evaluation_reuses_design=True)


def full_target_stratified_interval(data, alpha=.05, range_adaptive=False):
    """All-patient saturated plug-in with finite-sample binary-outcome coverage.

    Conditional on X,A, outcome error is one weighted sum of independent
    Bernoulli errors plus explicit half-count smoothing bias. Population design
    error is a mean of fixed true cell contrasts in [-1,1], bounded by Hoeffding.
    There is no nuisance fitting/validation split in this comparator.
    """
    if not 0 < alpha < 1:
        raise ValueError('Confidence failure must be in (0,1)')
    counts = data['target_training_counts'] + data['target_evaluation_counts']
    arm_counts = counts.sum(axis=3)
    cell_counts = arm_counts.sum(axis=2)
    total = cell_counts.sum(axis=1)
    masses = cell_counts / total[:, None]
    fitted_means = (counts[..., 1] + .5) / (arm_counts + 1)
    point = np.sum(masses * (fitted_means[:, :, 0] - fitted_means[:, :, 1]), axis=1)
    coefficients = np.repeat(masses[:, :, None], 2, axis=2).reshape(len(counts), -1)
    mean_failure = alpha / 2
    contrast_range = np.full(len(counts), 2.)
    if range_adaptive:
        range_failure = alpha / 10
        mean_failure = (alpha-range_failure)/2
        cells = counts.shape[1]
        lower, upper = binomial_bounds(counts[...,1],arm_counts,range_failure/(2*cells))
        contrast_lower = lower[:,:,0]-upper[:,:,1]
        contrast_upper = upper[:,:,0]-lower[:,:,1]
        contrast_range = np.minimum(2.,contrast_upper.max(axis=1)-contrast_lower.min(axis=1))
    radius = outcome_linear_radius(coefficients, arm_counts.reshape(len(counts), -1), mean_failure)
    radius += contrast_range * np.sqrt(np.log(2/mean_failure)/(2*total))
    interval = np.clip(np.column_stack((point - radius, point + radius)), -1, 1)
    return dict(interval=interval, point=point, contrast_range=contrast_range)

"""Conditional linear-balance scores with explicit model-error protection.

The base leave-out correction assumes linear conditional means at all sites.
Projection and mean-bias certificates extend the conditional bands to stated
approximation classes. This is separate from sparse/nonlinear production RoCE.
"""
import numpy as np


class UnsupportedDesign(ValueError):
    """The declared exact-balance model cannot be fit on these design rows."""


def linear_balance_design(features, target_mean, leverage_limit=.98):
    """Minimum-norm signed weights and leave-out variance geometry.

    features includes an intercept as its first column. target_mean includes
    the corresponding one. Bounds below use independent outcomes in [0,1].
    """
    features = np.asarray(features, dtype=float)
    target_mean = np.asarray(target_mean, dtype=float)
    if (features.ndim != 2 or features.shape[1] < 1 or target_mean.shape != (features.shape[1],)
            or not np.all(np.isfinite(features)) or not np.all(np.isfinite(target_mean))
            or not np.all(features[:, 0] == 1) or target_mean[0] != 1
            or not 0 < leverage_limit < 1):
        raise ValueError('Finite aligned features with an intercept are required')
    rows, columns = features.shape
    if rows <= columns:
        raise UnsupportedDesign('Exact balance requires more arm observations than columns')
    basis, triangular = np.linalg.qr(features, mode='reduced')
    if np.linalg.cond(triangular) > 1e8:
        raise UnsupportedDesign('Ill-conditioned balancing design')
    weights = basis.dot(np.linalg.solve(triangular.T, target_mean))
    leverage = np.sum(basis**2, axis=1)
    if leverage.max() >= leverage_limit:
        raise UnsupportedDesign('Leverage guard prevents unstable leave-out variance')
    balance_error = np.max(np.abs(features.T.dot(weights) - target_mean))
    if balance_error > 1e-9:
        raise ArithmeticError('Exact feature balance failed')
    diagonal = weights**2 / (1 - leverage)
    projected_diagonal = basis.T.dot(diagonal[:, None] * basis)
    correction_frobenius = (np.sum(diagonal**2 * (1 - 1.5 * leverage))
                            + .5 * np.sum(projected_diagonal**2))
    return dict(basis=basis, triangular=triangular, weights=weights, leverage=leverage, diagonal=diagonal,
                weight_square_sum=np.sum(weights**2), diagonal_square_sum=np.sum(diagonal**2),
                correction_frobenius=max(0., correction_frobenius),
                max_diagonal=diagonal.max(), balance_error=balance_error,
                negative_weight_fraction=np.mean(weights < 0))


def evaluate_linear_balance(design, outcomes):
    """Linear-balance score and leave-out variance (which may be negative).

    The variance is unbiased under linear conditional means; otherwise its
    projection bias needs the explicit certificate used by the caller.
    """
    outcomes = np.asarray(outcomes, dtype=float)
    if (outcomes.shape != design['weights'].shape or not np.all(np.isfinite(outcomes))
            or np.any((outcomes < 0) | (outcomes > 1))):
        raise ValueError('Aligned bounded outcomes are required')
    centered = outcomes - .5
    residual = centered - design['basis'].dot(design['basis'].T.dot(centered))
    variance = np.sum(design['diagonal'] * centered * residual)
    return dict(estimate=design['weights'].dot(outcomes), variance=variance,
                **{name: design[name] for name in ['weight_square_sum', 'diagonal_square_sum',
                  'correction_frobenius', 'max_diagonal', 'balance_error', 'negative_weight_fraction']})


def leaveout_dispersion_constants(summaries):
    """Base sub-gamma constants for exact linear conditional means."""
    size = summaries['estimate'].shape[-1]
    square_sum = summaries['weight_square_sum']
    diagonal = (size - 1) / size**2
    linear = 2 * np.max(square_sum, axis=-1) / size
    mean_frobenius = (np.sum(square_sum, axis=-1)**2
        + ((size - 1)**2 - 1) * np.sum(square_sum**2, axis=-1)) / size**4
    correction_frobenius = diagonal**2 * np.sum(summaries['correction_frobenius'], axis=-1)
    extra_linear = diagonal**2 * np.sum(summaries['diagonal_square_sum'], axis=-1) / 8
    constant = extra_linear + (mean_frobenius + correction_frobenius) / 8
    scale = diagonal * np.max(summaries['max_diagonal'], axis=-1) / 2
    return linear, constant, scale


def leaveout_dispersion_upper(summaries, alpha, variance_bias_sum=0., weighted_projection_error=0.):
    if not 0 < alpha < 1:
        raise ValueError('Failure allowance must be in (0,1)')
    if (np.any(np.asarray(variance_bias_sum) < 0) or np.any(np.asarray(weighted_projection_error) < 0)
            or not np.all(np.isfinite(variance_bias_sum)) or not np.all(np.isfinite(weighted_projection_error))):
        raise ValueError('Finite nonnegative approximation certificates required')
    estimates = summaries['estimate']
    size = estimates.shape[-1]
    if size == 1:
        return np.zeros(estimates.shape[:-1])
    observed = np.var(estimates, axis=-1) - (size - 1) / size**2 * np.sum(summaries['variance'], axis=-1)
    linear, constant, scale = leaveout_dispersion_constants(summaries)
    diagonal = (size - 1) / size**2
    squared_diagonal = np.sum(summaries['diagonal_square_sum'], axis=-1)
    constant = constant + diagonal**2/2 * (weighted_projection_error
        + np.sqrt(weighted_projection_error*squared_diagonal))
    logarithm = np.log(1 / alpha)
    shifted = observed + diagonal*variance_bias_sum + scale * logarithm
    discriminant = 2 * logarithm * linear * shifted + logarithm**2 * linear**2 + 2 * logarithm * constant
    return np.maximum(0., shifted + logarithm * linear + np.sqrt(np.maximum(0., discriminant)))


def contrast_summaries(arms):
    """Treatment groups are disjoint after conditioning on the full design."""
    result = {name: np.sum(value, axis=-1) for name, value in arms.items()
              if name not in ['estimate', 'max_diagonal', 'balance_error', 'negative_weight_fraction']}
    result['estimate'] = arms['estimate'][..., 0] - arms['estimate'][..., 1]
    result['max_diagonal'] = arms['max_diagonal'].max(axis=-1)
    return result


def conditional_linear_bands(source, target, valid_minimum, dispersion_mode='adaptive',
                             source_failure=.01, target_failure=.01, dispersion_failure=.02,
                             dispersion_budgets=None, mean_bias_bounds=None, target_bias_bound=0.):
    """Conditional TATE bands with optional valid class certificates.

    Mean-bias columns are (treated, control). Dispersion-budget keys are
    mu1, mu0 and contrast; their values are passed to leaveout_dispersion_upper.
    Certificate failure probabilities belong in the caller's total allocation.
    """
    if any(not 0 < value < 1 for value in [source_failure, target_failure, dispersion_failure]):
        raise ValueError('Each failure allowance must be in (0,1)')
    if not np.isfinite(target_bias_bound) or target_bias_bound < 0:
        raise ValueError('Target bias bound must be finite and nonnegative')
    if dispersion_mode not in ['arms', 'contrast', 'adaptive']:
        raise ValueError('Unknown dispersion mode')
    size = source['estimate'].shape[0]
    valid = np.broadcast_to(valid_minimum, (2,)).astype(float)
    if np.any(valid < 1) or np.any(valid > size) or np.any(valid != np.floor(valid)):
        raise ValueError('Integer valid-source guarantees in [1,K] required')
    intersection_count = int(valid.sum() - size)
    if dispersion_mode == 'contrast' and intersection_count < 1:
        raise ValueError('Contrast dispersion requires a guaranteed valid intersection')
    source_contrast = contrast_summaries(source)
    center = source_contrast['estimate'].mean()
    noise = np.sqrt(np.sum(source_contrast['weight_square_sum']) / (2 * size**2) * np.log(2/source_failure))
    target_point = target['estimate'][0] - target['estimate'][1]
    target_noise = np.sqrt(np.sum(target['weight_square_sum']) / 2 * np.log(2/target_failure))
    target_noise += target_bias_bound
    budgets = {} if dispersion_budgets is None else dispersion_budgets
    mean_bounds = np.zeros_like(source['estimate']) if mean_bias_bounds is None else np.asarray(mean_bias_bounds)
    if mean_bounds.shape != source['estimate'].shape or np.any(mean_bounds < 0) or not np.all(np.isfinite(mean_bounds)):
        raise ValueError('Nonnegative source-by-arm mean-bias bounds required')
    valid_bias = sum(np.sort(mean_bounds[:, arm])[-int(valid[arm]):].mean() for arm in range(2))
    arm_radius, contrast_radius = np.inf, np.inf
    use_arms = dispersion_mode != 'contrast'
    use_contrast = dispersion_mode != 'arms' and intersection_count > 0
    split = use_arms and use_contrast
    if use_arms:
        arm_radius = 0.
        for arm in range(2):
            if valid[arm] < size:
                summary = {name: value[:, arm] for name, value in source.items()}
                upper = leaveout_dispersion_upper(summary, dispersion_failure/(4 if split else 2),
                    **budgets.get(('mu1','mu0')[arm], {}))
                arm_radius += np.sqrt((size - valid[arm]) / valid[arm] * upper)
    if use_contrast:
        upper = leaveout_dispersion_upper(source_contrast, dispersion_failure/(2 if split else 1),
            **budgets.get('contrast', {}))
        contrast_radius = np.sqrt((size - intersection_count) / intersection_count * upper)
        contrast_radius += np.sort(mean_bounds.sum(axis=1))[-intersection_count:].mean()
    arm_radius += valid_bias
    bias = min(arm_radius, contrast_radius)
    source_band = np.array([center - noise - bias, center + noise + bias])
    target_band = np.array([target_point - target_noise, target_point + target_noise])
    combined = np.array([max(source_band[0], target_band[0]), min(source_band[1], target_band[1])])
    fallback = combined[0] > combined[1]
    if fallback:
        combined = source_band if np.ptp(source_band) < np.ptp(target_band) else target_band
    projection = np.clip(target_point, combined[0], combined[1])
    return dict(conditional_interval=combined, source_band=source_band, target_band=target_band,
                target_projection=projection, target_point=target_point, source_center=center,
                source_noise=noise, target_noise=target_noise, dispersion_radius=bias,
                arm_radius=arm_radius, contrast_radius=contrast_radius, fallback=fallback)


def continuous_balance_interval(source, target, valid_minimum, target_size,
                                dispersion_mode='adaptive', cate_range=2.):
    """Legacy v19 population interval with a declared CATE range bound."""
    if not 0 <= cate_range <= 2 or target_size < 1:
        raise ValueError('A declared CATE range in [0,2] and positive target size are required')
    result = conditional_linear_bands(source, target, valid_minimum, dispersion_mode)
    composition = cate_range*np.sqrt(np.log(200.)/(2*target_size))
    interval = np.clip(result['conditional_interval'] + composition*np.array([-1, 1]), -1, 1)
    baseline_noise = np.sqrt(np.sum(target['weight_square_sum']) / 2 * np.log(80.))
    baseline_composition = cate_range * np.sqrt(np.log(80.) / (2 * target_size))
    baseline = np.clip(result['target_point'] + (baseline_noise+baseline_composition)*np.array([-1, 1]), -1, 1)
    result.update(interval=interval, midpoint=interval.mean(),
                  target_projection=np.clip(result['target_projection'], interval[0], interval[1]),
                  target_interval=baseline, composition_radius=composition)
    return result

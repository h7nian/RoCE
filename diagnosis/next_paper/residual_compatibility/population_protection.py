"""Observed-data certificates for target population-composition uncertainty."""
import itertools
import numpy as np
from continuous_balance import linear_balance_design


def target_coefficient_certificate(features, treatment, outcomes, alpha=.005,
                                   support_lower=-1., support_upper=1., enumeration_limit=12):
    """One conditional sub-Gaussian ellipsoid protects range and sample variance.

    Support bounds must be declared independently of the observed sample.
    The target conditional outcome means must be linear in the supplied basis.
    """
    features = np.asarray(features, dtype=float)
    treatment, outcomes = np.asarray(treatment), np.asarray(outcomes, dtype=float)
    if features.ndim != 2 or features.shape[1] < 1 or len(features) < 2 or not 0 < alpha < 1:
        raise ValueError('A target feature matrix, n>=2 and valid alpha are required')
    if treatment.shape != (len(features),) or outcomes.shape != treatment.shape or not np.all(np.isin(treatment, [0, 1])):
        raise ValueError('Aligned binary treatment and outcome arrays required')
    if np.any((outcomes < 0) | (outcomes > 1)) or not np.all(np.isfinite(outcomes)):
        raise ValueError('Outcomes must be bounded in [0,1]')
    lower = np.broadcast_to(support_lower, (features.shape[1],))
    upper = np.broadcast_to(support_upper, (features.shape[1],))
    if not np.all(np.isfinite(features)) or not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)) or np.any(lower >= upper) or np.any(features < lower) or np.any(features > upper):
        raise ValueError('Observed features must lie within finite declared support bounds')
    matrix = np.column_stack((np.ones(len(features)), features))
    target_mean = matrix.mean(axis=0)
    coefficients, proxies = [], []
    for arm in [1, 0]:
        selected = treatment == arm
        design = linear_balance_design(matrix[selected], target_mean)
        inverse = np.linalg.solve(design['triangular'], np.eye(matrix.shape[1]))
        coefficients.append(inverse.dot(design['basis'].T.dot(outcomes[selected])))
        proxies.append(.25*inverse.dot(inverse.T))
    estimate = (coefficients[0]-coefficients[1])[1:]
    proxy = (proxies[0]+proxies[1])[1:, 1:]
    dimension = len(estimate)
    logarithm = np.log(1/alpha)
    radius = np.sqrt(dimension+2*np.sqrt(dimension*logarithm)+2*logarithm)
    widths = upper-lower
    if dimension <= enumeration_limit:
        directions = np.array(list(itertools.product([-1., 1.], repeat=dimension)))*widths
        range_upper = np.max(directions.dot(estimate)+radius*np.sqrt(
            np.einsum('ij,jk,ik->i', directions, proxy, directions)))
    else:
        # A safe upper bound avoids exponential enumeration at larger p.
        range_upper = np.sum(widths*(np.abs(estimate)+radius*np.sqrt(np.diag(proxy))))
    range_upper = min(2., max(0., float(range_upper)))
    feature_covariance = np.atleast_2d(np.cov(features, rowvar=False, ddof=1))
    factor = np.linalg.cholesky(proxy)
    uncertainty = np.linalg.eigvalsh(factor.T.dot(feature_covariance).dot(factor)).max()
    estimated_variance = max(0., estimate.dot(feature_covariance).dot(estimate))
    variance_upper = (np.sqrt(estimated_variance)+radius*np.sqrt(max(0., uncertainty)))**2
    variance_upper = min(variance_upper, len(features)/(len(features)-1)*range_upper**2/4)
    return dict(coefficient=estimate, proxy=proxy, ellipsoid_radius=radius,
                range_upper=range_upper, sample_variance_upper=variance_upper,
                feature_covariance=feature_covariance, failure=alpha)


def composition_radius(certificate, sample_size, alpha, method='bernstein'):
    if sample_size < 2 or not 0 < alpha < 1:
        raise ValueError('n>=2 and a failure allowance in (0,1) are required')
    if method == 'universal':
        return 2*np.sqrt(np.log(2/alpha)/(2*sample_size))
    if method == 'range':
        return certificate['range_upper']*np.sqrt(np.log(2/alpha)/(2*sample_size))
    if method == 'bernstein':
        logarithm = np.log(4/alpha)
        return (np.sqrt(2*certificate['sample_variance_upper']*logarithm/sample_size)
                + 7*certificate['range_upper']*logarithm/(3*(sample_size-1)))
    if method == 'adaptive':
        return min(composition_radius(certificate, sample_size, alpha/2, 'range'),
                   composition_radius(certificate, sample_size, alpha/2, 'bernstein'))
    raise ValueError('Unknown population-composition method')


def expand_conditional_interval(interval, radius):
    interval = np.asarray(interval)
    if interval.shape[-1] != 2 or np.any(interval[..., 0] > interval[..., 1]) or np.any(np.asarray(radius) < 0):
        raise ValueError('Ordered intervals and nonnegative radii required')
    return np.clip(interval + np.asarray(radius)[..., None]*np.array([-1., 1.]), -1, 1)


def projection_error_budget(summary, residual_l2_bounds, max_nonlinear=None):
    """Worst unknown-subset certificate for leave-out correction misspecification."""
    bounds = np.broadcast_to(residual_l2_bounds, summary['estimate'].shape)
    if np.any(bounds < 0) or not np.all(np.isfinite(bounds)):
        raise ValueError('Finite nonnegative residual norm bounds required')
    size = bounds.shape[-1]
    count = size if max_nonlinear is None else max_nonlinear
    if not 0 <= count <= size or count != int(count):
        raise ValueError('An integer upper bound on nonlinear source count is required')
    def largest_sum(values):
        return np.sum(np.sort(values, axis=-1)[..., -int(count):], axis=-1) if count else np.zeros(values.shape[:-1])
    centering = .5*np.sqrt(summary['diagonal_square_sum'])*bounds
    if 'projection_commutator_norm' in summary:
        refined = (summary['max_diagonal']*bounds**2 + .5*np.sqrt(summary['sample_size'])
                   * summary['projection_commutator_norm']*bounds)
        centering = np.minimum(centering, refined)
    return dict(variance_bias_sum=largest_sum(centering),
                weighted_projection_error=largest_sum(summary['max_diagonal']**2*bounds**2))


def projection_geometry(design):
    """Norm of (I-H) diag(d) H, computed through the thin QR basis."""
    basis, diagonal = design['basis'], design['diagonal']
    weighted = basis.T.dot(diagonal[:, None]*basis)
    square = basis.T.dot(diagonal[:, None]**2*basis)-weighted.dot(weighted)
    square = (square+square.T)/2
    return dict(projection_commutator_norm=np.sqrt(max(0., np.linalg.eigvalsh(square).max())),
                sample_size=len(diagonal), weight_absolute_sum=np.sum(np.abs(design['weights'])))

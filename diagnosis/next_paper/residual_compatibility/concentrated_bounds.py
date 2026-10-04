"""Direct bounds for pooled moments and finite-stratum calibration remainders.

The quadratic bound uses patient-level Gaussian MGF domination, not a Gaussian
approximation to the fitted score means. See sections/hybrid.tex for the proof.
"""
import math
import numpy as np


def _moment_inputs(means, mean_variances, sample_sizes, alpha):
    means = np.asarray(means, dtype=float)
    if means.ndim != 2 or min(means.shape) < 1 or not 0 < alpha < 1:
        raise ValueError("Nonempty repeat-by-source matrix and alpha in (0,1) required")
    variances = np.broadcast_to(mean_variances, means.shape)
    sizes = np.broadcast_to(sample_sizes, means.shape)
    if any(not np.all(np.isfinite(value)) for value in [means, variances, sizes]):
        raise ValueError("Finite means, variances, and sample sizes required")
    if np.any(variances < -1e-12) or np.any(sizes < 2) or np.any(sizes != np.floor(sizes)):
        raise ValueError("Unbiased nonnegative mean variances and integer sizes >=2 required")
    return means, np.maximum(variances, 0), sizes


def pooled_mean_interval(means, mean_variances, sample_sizes, absolute_bounds, alpha):
    """Empirical Bernstein for independent, possibly nonidentical patients.

    Equal source weights require equal site sizes here. mean_variances contains
    unbiased S_j^2/n_j computed from the same observations as means. The supplied
    bounds limit absolute score values, not merely their within-site ranges.
    """
    means, variances, sizes = _moment_inputs(means, mean_variances, sample_sizes, alpha)
    if np.any(sizes != sizes[:, :1]):
        raise ValueError("Pooled patient calibration requires equal site sizes")
    bounds = np.broadcast_to(absolute_bounds, means.shape)
    if not np.all(np.isfinite(bounds)) or np.any(bounds < 0):
        raise ValueError("Finite nonnegative absolute score bounds required")
    total = means.shape[1] * sizes[:, 0]
    center = means.mean(axis=1)
    within = np.sum(sizes * (sizes - 1) * variances, axis=1)
    between = np.sum(sizes * (means - center[:, None])**2, axis=1)
    pooled_variance = (within + between) / (total * (total - 1))
    score_range = 2 * np.max(bounds, axis=1)
    logarithm = math.log(4 / alpha)
    radius = (np.sqrt(2 * pooled_variance * logarithm)
              + 7 * score_range * logarithm / (3 * (total - 1)))
    return np.column_stack((center - radius, center + radius))


def exponential_dispersion_constants(score_ranges, sample_sizes, shape):
    ranges = np.broadcast_to(score_ranges, shape)
    sizes = np.broadcast_to(sample_sizes, shape)
    if (np.any(sizes < 2) or np.any(sizes != np.floor(sizes)) or np.any(ranges < 0)
            or not np.all(np.isfinite(ranges)) or not np.all(np.isfinite(sizes))):
        raise ValueError("Finite ranges and integer sizes >=2 required")
    count = shape[1]
    proxy = ranges**2 / (4 * sizes)
    linear = 4 * np.max(proxy, axis=1) / count
    frobenius = (np.sum(proxy, axis=1)**2 + np.sum(
        (sizes / (sizes - 1) * (count - 1)**2 - 1) * proxy**2, axis=1)) / count**4
    negative_eigenvalue = np.max((count - 1) / count**2 * proxy / (sizes - 1), axis=1)
    return linear, 2 * frobenius, 2 * negative_eigenvalue


def exponential_dispersion_upper(means, mean_variances, sample_sizes, score_ranges, alpha):
    means, variances, sizes = _moment_inputs(means, mean_variances, sample_sizes, alpha)
    count = means.shape[1]
    linear, quadratic, scale = exponential_dispersion_constants(score_ranges, sizes, means.shape)
    if count == 1:
        return np.zeros(means.shape[0])
    observed = (np.mean((means - means.mean(axis=1)[:, None])**2, axis=1)
                - (count - 1) / count**2 * np.sum(variances, axis=1))
    logarithm = math.log(1 / alpha)
    shifted = observed + scale * logarithm
    discriminant = 2 * logarithm * linear * shifted + logarithm**2 * linear**2 + 2 * logarithm * quadratic
    return np.maximum(0, shifted + logarithm * linear + np.sqrt(np.maximum(0, discriminant)))


def concentrated_dispersion_interval(means, mean_variances, sample_sizes, score_ranges,
                                     absolute_bounds, valid_minimum, alpha, average_bias=0.):
    means = np.asarray(means)
    count = means.shape[1]
    if not 1 <= valid_minimum <= count or valid_minimum != math.floor(valid_minimum):
        raise ValueError("An integer guaranteed source count in [1,K] is required")
    mean_alpha = alpha if valid_minimum == count else alpha / 2
    interval = pooled_mean_interval(means, mean_variances, sample_sizes, absolute_bounds, mean_alpha)
    radius = np.broadcast_to(average_bias, (means.shape[0],)).copy()
    if np.any(radius < 0) or not np.all(np.isfinite(radius)):
        raise ValueError("Finite nonnegative average-drift bounds required")
    if valid_minimum < count:
        energy = exponential_dispersion_upper(means, mean_variances, sample_sizes, score_ranges, alpha / 2)
        radius += np.sqrt((count - valid_minimum) / valid_minimum * energy)
    interval[:, 0] -= radius
    interval[:, 1] += radius
    return interval


def simplex_linear_max(coefficients,lower,upper):
    """Exact LP over box-constrained probabilities with total mass one."""
    coefficients,lower,upper=np.broadcast_arrays(coefficients,lower,upper)
    order=np.argsort(-coefficients,axis=-1)
    mass=lower.copy();remaining=1-np.sum(mass,axis=-1)
    if np.any(remaining<-1e-9) or np.any(np.sum(upper,axis=-1)<1-1e-9):
        raise ValueError("Empty probability region")
    for rank in range(coefficients.shape[-1]):
        index=order[...,rank:rank+1]
        capacity=np.take_along_axis(upper-mass,index,axis=-1)[...,0]
        addition=np.minimum(np.maximum(remaining,0),np.maximum(capacity,0))
        np.put_along_axis(mass,index,np.take_along_axis(mass,index,axis=-1)+addition[...,None],axis=-1)
        remaining-=addition
    if np.any(remaining>1e-8):raise ArithmeticError("Probability LP did not allocate unit mass")
    return np.sum(coefficients*mass,axis=-1)


def grouped_probability_max(coefficients,lower,upper,group_lower,group_upper):
    """Exact linear maximum for (cell,arm) probabilities and cell-mass bounds."""
    coefficients, lower, upper = np.broadcast_arrays(coefficients, lower, upper)
    if coefficients.shape[-1] != 2:
        raise ValueError("Grouped probability solver expects two arms per cell")
    mass = np.array(lower, dtype=float, copy=True)
    groups = mass.shape[-2]
    group_min=np.maximum(group_lower,lower.sum(axis=-1))
    group_max=np.minimum(group_upper,upper.sum(axis=-1))
    for group in range(groups):
        remaining=np.maximum(0,group_min[...,group]-mass[...,group,:].sum(axis=-1))
        order=np.argsort(-coefficients[...,group,:],axis=-1)
        for rank in [0,1]:
            index=order[...,rank:rank+1]
            block=mass[...,group,:]
            capacity=np.take_along_axis(upper[...,group,:]-block,index,axis=-1)[...,0]
            addition=np.minimum(remaining,np.maximum(capacity,0))
            np.put_along_axis(block,index,np.take_along_axis(block,index,axis=-1)+addition[...,None],axis=-1)
            remaining-=addition
    flat=mass.reshape(mass.shape[:-2]+(-1,))
    flat_upper=upper.reshape(upper.shape[:-2]+(-1,))
    flat_coefficients=coefficients.reshape(coefficients.shape[:-2]+(-1,))
    remaining=1-flat.sum(axis=-1)
    order=np.argsort(-flat_coefficients,axis=-1)
    for rank in range(2*groups):
        index=order[...,rank:rank+1]
        group=index[...,0]//2
        group_capacity=np.take_along_axis(group_max-mass.sum(axis=-1),group[...,None],axis=-1)[...,0]
        capacity=np.take_along_axis(flat_upper-flat,index,axis=-1)[...,0]
        addition=np.minimum(np.maximum(remaining,0),np.maximum(0,np.minimum(capacity,group_capacity)))
        np.put_along_axis(flat,index,np.take_along_axis(flat,index,axis=-1)+addition[...,None],axis=-1)
        remaining-=addition
    if np.any(abs(remaining)>1e-8):raise ValueError("Inconsistent grouped probability bounds")
    return np.sum(flat_coefficients*flat,axis=-1)

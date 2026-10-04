"""Joint TATE score noise and conditional nuisance certificates.

Each group contains independent held-out patients. Centers, score supports and
weights must be fixed before those outcomes/designs are evaluated. The source
valid sets may differ by arm; only the noise and nuisance contrast are joined.
"""
import numpy as np
from conditional_remainders import (conditional_design_coefficients, outcome_linear_radius,
                                    full_target_stratified_interval)
from concentrated_bounds import exponential_dispersion_upper


def weighted_empirical_mean_radius(means, mean_variances, sizes, weights,
                                   lower, upper, centers, alpha):
    """Maurer--Pontil Theorem 11 after scaling independent patient observations.

    Means/mean_variances are group means and unbiased S_g^2/n_g. Scaling each
    centered patient by N*w_g/n_g makes its overall mean the desired weighted
    sum. The pooled empirical variance includes nonidentical-group means.
    """
    means = np.asarray(means, dtype=float)
    if means.ndim != 2 or min(means.shape) == 0 or not 0 < alpha < 1:
        raise ValueError('A nonempty repeat-by-group matrix and alpha in (0,1) required')
    variances, sizes, weights, lower, upper, centers = [np.broadcast_to(x, means.shape)
        for x in [mean_variances, sizes, weights, lower, upper, centers]]
    if any(not np.all(np.isfinite(x)) for x in [means,variances,sizes,weights,lower,upper,centers]):
        raise ValueError('Finite group statistics and predictable bounds required')
    if np.any(sizes < 2) or np.any(sizes != np.floor(sizes)) or np.any(variances < -1e-12) or np.any(lower > upper):
        raise ValueError('Integer sizes >=2, nonnegative variances, ordered bounds required')
    total = sizes.sum(axis=1)
    centered = means - centers
    weighted_sum = np.sum(weights * centered, axis=1)
    squares = np.sum(weights**2 * ((sizes-1)/sizes*np.maximum(variances,0) + centered**2/sizes), axis=1)
    empirical_variance = np.maximum(0,(total*squares-weighted_sum**2)/(total-1))
    scaled_lower = weights*(lower-centers)/sizes
    scaled_upper = weights*(upper-centers)/sizes
    score_range = np.max(np.maximum(scaled_lower,scaled_upper),axis=1) - np.min(np.minimum(scaled_lower,scaled_upper),axis=1)
    logarithm = np.log(4/alpha)
    radius = np.sqrt(2*empirical_variance*logarithm) + 7*total*score_range*logarithm/(3*(total-1))
    return radius, empirical_variance


def noise_radius(data, alpha, mode='joint', centering='training'):
    if mode not in ['joint','separate'] or centering not in ['training','zero']:
        raise ValueError('Unknown contrast-noise mode or predictable centering')
    repeats,count = data['source'].shape[:2]
    target_center = data['target_training_prediction_mean'] if centering=='training' else np.zeros(repeats)
    source_center = data['source_training_mean'] if centering=='training' else np.zeros_like(data['source'])
    target_bounds = data['target_prediction_bounds']
    target_arguments = [data['target'][:,0,None],data['target_sample_covariance'][:,0,0,None],
        data['target_sample_size'],1.,target_bounds[:,0,None],target_bounds[:,1,None],target_center[:,None]]
    if mode=='separate':
        radius = weighted_empirical_mean_radius(*target_arguments,alpha/3)[0]
        for arm in [0,1]:
            radius += weighted_empirical_mean_radius(data['source'][:,:,arm],data['sample_variance'][:,:,arm],
                data['sample_sizes'],1/count,data['source_score_lower'][:,:,arm],data['source_score_upper'][:,:,arm],
                source_center[:,:,arm],alpha/3)[0]
        return radius
    contrast = data['source'][:,:,0]-data['source'][:,:,1]
    variance = data['sample_variance'].sum(axis=2)-2*data['source_sample_cross_covariance']
    sizes = np.broadcast_to(data['sample_sizes'],contrast.shape)
    means = np.column_stack((data['target'][:,0],contrast))
    variances = np.column_stack((data['target_sample_covariance'][:,0,0],variance))
    group_sizes = np.column_stack((np.full(repeats,data['target_sample_size']),sizes))
    weights = np.r_[1.,np.full(count,1/count)]
    lower = np.column_stack((target_bounds[:,0],data['source_contrast_lower']))
    upper = np.column_stack((target_bounds[:,1],data['source_contrast_upper']))
    centers = np.column_stack((target_center,source_center[:,:,0]-source_center[:,:,1]))
    return weighted_empirical_mean_radius(means,variances,group_sizes,weights,lower,upper,centers,alpha)[0]


def nuisance_radius(data, coefficients, outcome_failure, mode='joint'):
    bounds = coefficients['source_coefficient_bound']
    counts = data['target_training_counts'].sum(axis=3)
    if mode=='joint':
        return outcome_linear_radius(bounds.reshape(len(bounds),-1),counts.reshape(len(counts),-1),outcome_failure)
    if mode=='separate':
        return sum(outcome_linear_radius(bounds[:,:,arm],counts[:,:,arm],outcome_failure/2) for arm in [0,1])
    raise ValueError('Unknown nuisance contrast mode')


def infer_joint_contrast(data, valid_minimum, noise_mode='joint', nuisance_mode='joint',
                         source_calibration='empirical_design', envelope=None, coefficients=None,
                         noise_failure=.025, dispersion_failure=.01, nuisance_failure=.005,
                         target_failure=.01, centering='training', dispersion_mode='arms'):
    if data.get('nuisance_structure')!='saturated_design_weights':
        raise ValueError('Joint conditional inference requires the declared design-only saturated weights')
    allowances = np.array([noise_failure,dispersion_failure,nuisance_failure,target_failure])
    if np.any(allowances<=0) or not np.all(np.isfinite(allowances)) or allowances.sum()>.05+1e-14:
        raise ValueError('Positive confidence allowances must sum to at most .05')
    repeats,count = data['source'].shape[:2]
    valid = np.broadcast_to(valid_minimum,(2,))
    if np.any(valid<1) or np.any(valid>count) or np.any(valid!=np.floor(valid)):
        raise ValueError('Integer valid-source guarantees required in each arm')
    if dispersion_mode not in ['arms','contrast']:
        raise ValueError('Unknown source dispersion mode')
    certificate_valid = valid
    if dispersion_mode == 'contrast':
        overlap_guarantee = int(valid.sum()-count)
        if overlap_guarantee < 1:
            raise ValueError('Contrast dispersion requires g1+g0>K; keep arm dispersion otherwise')
        certificate_valid = np.full(2,overlap_guarantee)
    if coefficients is None:
        coefficients = conditional_design_coefficients(data,certificate_valid,nuisance_failure*.2,nuisance_failure*.4,
                                                        source_calibration,envelope)
    if (not np.isclose(coefficients['target_design_failure'],nuisance_failure*.2,atol=1e-15,rtol=1e-12)
            or not np.isclose(coefficients['source_design_failure'],nuisance_failure*.4,atol=1e-15,rtol=1e-12)
            or np.any(coefficients['valid_minimum']!=certificate_valid)
            or coefficients['source_calibration']!=source_calibration):
        raise ValueError('Precomputed coefficient certificate has a different confidence contract')
    bias_radius = nuisance_radius(data,coefficients,nuisance_failure*.4,nuisance_mode)
    score_radius = noise_radius(data,noise_failure,noise_mode,centering)
    dispersion_radius = np.zeros(repeats)
    if dispersion_mode == 'arms':
        energies = np.zeros((repeats,2))
        for arm in [0,1]:
            if valid[arm]<count:
                energies[:,arm] = exponential_dispersion_upper(data['source'][:,:,arm],data['sample_variance'][:,:,arm],
                    data['sample_sizes'],data['score_ranges'][:,:,arm],dispersion_failure/2)
                dispersion_radius += np.sqrt((count-valid[arm])/valid[arm]*energies[:,arm])
    else:
        contrast = data['source'][:,:,0]-data['source'][:,:,1]
        variance = data['sample_variance'].sum(axis=2)-2*data['source_sample_cross_covariance']
        energy = exponential_dispersion_upper(contrast,variance,data['sample_sizes'],
            data['source_contrast_upper']-data['source_contrast_lower'],dispersion_failure)
        energies = energy[:,None]
        dispersion_radius = np.sqrt((count-overlap_guarantee)/overlap_guarantee*energy)
    center = data['target'][:,0]+np.mean(data['source'][:,:,0]-data['source'][:,:,1],axis=1)
    radius = score_radius+bias_radius+dispersion_radius
    source_interval = np.clip(np.column_stack((center-radius,center+radius)),-1,1)
    target = full_target_stratified_interval(data,alpha=target_failure)
    target_interval = target['interval']
    intersection = np.column_stack((np.maximum(source_interval[:,0],target_interval[:,0]),
                                    np.minimum(source_interval[:,1],target_interval[:,1])))
    empty = intersection[:,0]>intersection[:,1]
    source_shorter = np.diff(source_interval,axis=1)[:,0] < np.diff(target_interval,axis=1)[:,0]
    fallback = np.where(empty,np.where(source_shorter,2,1),0)
    fallback_interval = np.where(source_shorter[:,None],source_interval,target_interval)
    interval = np.where(empty[:,None],fallback_interval,intersection)
    return dict(interval=interval,source_interval=source_interval,target_interval=target_interval,
        point_midpoint=interval.mean(axis=1),point_source=np.clip(center,interval[:,0],interval[:,1]),
        source_center=center,noise_radius=score_radius,nuisance_radius=bias_radius,
        dispersion_radius=dispersion_radius,dispersion_upper=energies,fallback=fallback,
        dispersion_mode=dispersion_mode,certificate_valid_minimum=certificate_valid,
        length_envelope=np.minimum(np.diff(source_interval,axis=1)[:,0],np.diff(target_interval,axis=1)[:,0]))

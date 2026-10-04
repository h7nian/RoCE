"""Exact design calibration benchmark for finite binary-outcome strata.

This changes the score construction and sampling argument. It is NOT a drop-in
high-dimensional RoCE fit: positive cell counts and saturated strata are used.
"""
import numpy as np
from scipy.stats import norm
from fitted_strata import binomial_bounds


def block_dispersion_constants(coefficients, counts):
    """Gaussian MGF proxy constants for sums of independent cell means.

    Each cell contains Bernoulli outcomes with Hoeffding variance proxy 1/4.
    The statistic subtracts unbiased within-cell variance estimates, giving a
    zero-diagonal patient quadratic form even with different cell expectations.
    """
    coefficients, counts = np.broadcast_arrays(coefficients, counts)
    if coefficients.ndim!=3 or np.any(counts<2) or np.any(counts!=np.floor(counts)):
        raise ValueError('Repeat-by-site-by-cell arrays with integer counts >=2 required')
    if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(counts)):
        raise ValueError('Finite coefficients and counts required')
    sites=counts.shape[1]
    block_proxy=coefficients**2/(4*counts)
    site_proxy=block_proxy.sum(axis=2)
    linear=4*site_proxy.max(axis=1)/sites
    constant_part=site_proxy.sum(axis=1)**2+((sites-1)**2-1)*np.sum(site_proxy**2,axis=1)
    negative_part=(sites-1)**2*np.sum(block_proxy**2/(counts-1),axis=(1,2))
    quadratic=2*(constant_part+negative_part)/sites**4
    scale=2*(sites-1)/sites**2*np.max(block_proxy/(counts-1),axis=(1,2))
    return linear,quadratic,scale,site_proxy


def block_dispersion_upper(means, mean_variances, coefficients, counts, alpha):
    if not 0<alpha<1:
        raise ValueError('Failure allowance must be in (0,1)')
    linear,quadratic,scale,_=block_dispersion_constants(coefficients,counts)
    sites=means.shape[1]
    if sites==1:return np.zeros(len(means))
    observed=np.var(means,axis=1)-(sites-1)/sites**2*np.sum(mean_variances,axis=1)
    logarithm=np.log(1/alpha)
    shifted=observed+scale*logarithm
    return np.maximum(0,shifted+logarithm*linear+np.sqrt(np.maximum(0,
        2*logarithm*linear*shifted+logarithm**2*linear**2+2*logarithm*quadratic)))


def target_stratified_components(target_counts, range_failure, outcome_failure, composition_failure):
    """Unsmoothed cell means when observed; explicit bias protects empty arms."""
    if any(not 0<failure<1 for failure in [range_failure,outcome_failure,composition_failure]):
        raise ValueError('All confidence allowances must be in (0,1)')
    if target_counts.ndim!=4 or target_counts.shape[2:]!=(2,2) or np.any(target_counts<0):
        raise ValueError('Nonnegative binary-outcome cell counts required')
    counts=target_counts.sum(axis=3)
    total=counts.sum(axis=(1,2))
    if np.any(total<1):raise ValueError('Target samples must be nonempty')
    masses=counts.sum(axis=2)/total[:,None]
    observed=counts>0
    means=np.where(observed,target_counts[...,1]/np.maximum(counts,1),.5)
    point=np.sum(masses*(means[:,:,0]-means[:,:,1]),axis=1)
    empty_bias=.5*np.sum(masses[:,:,None]*(~observed),axis=(1,2))
    proxy=np.sum(masses[:,:,None]**2*observed/np.maximum(counts,1),axis=(1,2))/4
    outcome_radius=empty_bias+np.sqrt(2*proxy*np.log(2/outcome_failure))
    lower,upper=binomial_bounds(target_counts[...,1],counts,range_failure/(2*counts.shape[1]))
    contrast_range=np.minimum(2.,(upper[:,:,0]-lower[:,:,1]).max(axis=1)-(lower[:,:,0]-upper[:,:,1]).min(axis=1))
    composition_radius=contrast_range*np.sqrt(np.log(2/composition_failure)/(2*total))
    return dict(point=point,masses=masses,outcome_radius=outcome_radius,
                composition_radius=composition_radius,contrast_range=contrast_range)


def conditional_standardization_interval(data, valid_minimum, cohort='all', dispersion_mode='adaptive', point_rule='midpoint'):
    if cohort not in ['all','heldout'] or dispersion_mode not in ['arms','contrast','adaptive'] or point_rule not in ['midpoint','target_projection']:
        raise ValueError('Unknown cohort or dispersion certificate')
    source_counts=data['source_held_counts'].astype(np.int32)
    if cohort=='all':source_counts=source_counts+data['source_training_counts']
    sizes=source_counts.sum(axis=4)
    repeats,sites,cells,arms=sizes.shape
    valid=np.broadcast_to(valid_minimum,(2,)).astype(float)
    if arms!=2 or np.any(valid<1) or np.any(valid>sites) or np.any(valid!=np.floor(valid)):
        raise ValueError('Two arms and integer valid-source guarantees required')
    overlap=int(valid.sum()-sites)
    if dispersion_mode=='contrast' and overlap<1:
        raise ValueError('Contrast certificate requires a guaranteed valid intersection')
    target_counts=data['target_training_counts']+data['target_evaluation_counts']
    # One range event, one population-composition event, and conditional source,
    # target and dispersion events. Total failure allowance is exactly .05.
    range_failure,composition_failure=.005,.01
    source_failure,target_failure,dispersion_failure=.01,.01,.015
    target=target_stratified_components(target_counts,range_failure,target_failure,composition_failure)
    baseline=target_stratified_components(target_counts,.005,.0225,.0225)
    target_interval=np.column_stack((target['point']-target['outcome_radius'],target['point']+target['outcome_radius']))
    available=np.all(sizes>=2,axis=(1,2,3))
    safe_sizes=np.maximum(sizes,2)
    fractions=source_counts[...,1]/safe_sizes
    masses=target['masses']
    coefficients=np.broadcast_to(masses[:,None,:,None],sizes.shape)
    estimates=np.sum(coefficients*fractions,axis=2)
    mean_variance=np.sum(coefficients**2*fractions*(1-fractions)/(safe_sizes-1),axis=2)
    arm_proxy=np.sum(coefficients**2/(4*safe_sizes),axis=2)
    center=np.mean(estimates[:,:,0]-estimates[:,:,1],axis=1)
    noise_radius=np.sqrt(2*np.sum(arm_proxy,axis=(1,2))/sites**2*np.log(2/source_failure))
    arm_radius=np.full(repeats,np.inf)
    contrast_radius=np.full(repeats,np.inf)
    contrast_available=overlap>=1 and dispersion_mode!='arms'
    arm_available=dispersion_mode!='contrast'
    split=arm_available and contrast_available
    if arm_available:
        arm_radius=np.zeros(repeats)
        alpha=dispersion_failure/(4 if split else 2)
        for arm in [0,1]:
            if valid[arm]<sites:
                energy=block_dispersion_upper(estimates[:,:,arm],mean_variance[:,:,arm],
                    coefficients[:,:,:,arm],safe_sizes[:,:,:,arm],alpha)
                arm_radius+=np.sqrt((sites-valid[arm])/valid[arm]*energy)
    if contrast_available:
        signed=coefficients.copy();signed[:,:,:,1]*=-1
        energy=block_dispersion_upper(estimates[:,:,0]-estimates[:,:,1],mean_variance.sum(axis=2),
            signed.reshape(repeats,sites,-1),safe_sizes.reshape(repeats,sites,-1),
            dispersion_failure/(2 if split else 1))
        contrast_radius=np.sqrt((sites-overlap)/overlap*energy)
    bias_radius=np.minimum(arm_radius,contrast_radius)
    source_interval=np.column_stack((center-noise_radius-bias_radius,center+noise_radius+bias_radius))
    intersection=np.column_stack((np.maximum(source_interval[:,0],target_interval[:,0]),
                                   np.minimum(source_interval[:,1],target_interval[:,1])))
    empty=intersection[:,0]>intersection[:,1]
    source_shorter=np.diff(source_interval,axis=1)[:,0]<np.diff(target_interval,axis=1)[:,0]
    fallback_interval=np.where(source_shorter[:,None],source_interval,target_interval)
    cohort_interval=np.where(empty[:,None],fallback_interval,intersection)
    cohort_interval=np.where(available[:,None],cohort_interval,target_interval)
    interval=np.clip(cohort_interval+target['composition_radius'][:,None]*np.array([-1.,1.]),-1,1)
    baseline_radius=baseline['outcome_radius']+baseline['composition_radius']
    baseline_interval=np.clip(np.column_stack((baseline['point']-baseline_radius,baseline['point']+baseline_radius)),-1,1)
    # Feasible but unprotected comparison: ignores source transport bias.
    source_variance=np.sum(mean_variance,axis=(1,2))/sites**2
    cell_contrast=fractions[:,:,:,0].mean(axis=1)-fractions[:,:,:,1].mean(axis=1)
    target_size=target_counts.sum(axis=(1,2,3))
    composition_variance=np.sum(masses*(cell_contrast-center[:,None])**2,axis=1)/np.maximum(target_size-1,1)
    naive_radius=norm.isf(.025)*np.sqrt(np.maximum(0,source_variance+composition_variance))
    naive_interval=np.clip(np.column_stack((center-naive_radius,center+naive_radius)),-1,1)
    naive_interval=np.where((available&(target_size>=2))[:,None],naive_interval,baseline_interval)
    midpoint=interval.mean(axis=1)
    projection=np.clip(target['point'],cohort_interval[:,0],cohort_interval[:,1])
    projection=np.clip(projection,interval[:,0],interval[:,1])
    return dict(interval=interval,point_midpoint=midpoint,point_target_projection=projection,
        point=midpoint if point_rule=='midpoint' else projection,point_rule=point_rule,
        cohort_interval=cohort_interval,source_cohort_interval=np.where(available[:,None],source_interval,np.nan),
        source_center=np.where(available,center,np.nan),
        target_point=baseline['point'],target_interval=baseline_interval,naive_interval=naive_interval,noise_radius=np.where(available,noise_radius,np.nan),
        dispersion_radius=np.where(available,bias_radius,np.nan),composition_radius=target['composition_radius'],
        arm_dispersion_radius=arm_radius,contrast_dispersion_radius=contrast_radius,
        unavailable=~available,fallback=empty&available,source_arm_estimates=np.where(available[:,None,None],estimates,np.nan))

"""Exact Bernoulli moment matching for local weak-bias interval lower bounds."""
import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.stats import binom


def mixture_chi_square(sample_size, shift, valid_probability=.8, information_factor=4.):
    """Chi-square of a mean-matched Bernoulli site mixture against the null.

    information_factor=4 uses N treated observations; factor=2 uses N i.i.d.
    (A,Y) observations with independent A~Bernoulli(1/2).

    The first-order likelihood term cancels. A positive series from degree two
    avoids subtracting nearly equal numbers in the three-term closed formula.
    """
    if (sample_size < 2 or sample_size != int(sample_size) or not 0 < valid_probability < 1
            or information_factor not in [2.,4.]):
        raise ValueError('An integer sample size >=2 and a mixture probability required')
    ratio = valid_probability/(1-valid_probability)
    if not np.isfinite(shift) or shift < 0 or max(shift, ratio*shift) >= .5:
        raise ValueError('Both Bernoulli probabilities must be interior')
    if shift == 0:
        return 0.
    degree = np.arange(2, int(sample_size)+1)
    log_moment, _ = logsumexp(np.vstack((np.full(len(degree), np.log(valid_probability)),
        np.log1p(-valid_probability)+degree*np.log(ratio))),
        b=np.vstack((np.ones(len(degree)), (-1.)**degree)), axis=0, return_sign=True)
    terms = (gammaln(sample_size+1)-gammaln(degree+1)-gammaln(sample_size-degree+1)
             +degree*np.log(information_factor*shift**2)+2*log_moment)
    return float(np.exp(logsumexp(terms)))


def local_length_lower_bound(sample_size, target_size, source_count, shift_constant=.1,
                             valid_probability=.8, valid_fraction=.75, alpha=.05,
                             experiment='fixed_arm'):
    """All-valid expected-length lower bound for honest local-class intervals.

    The alternative has linear constant site means, bounded per-source weak
    shifts and at least ceil(valid_fraction*K) valid sites after conditioning.
    Sizes count treated patients in fixed_arm mode and all patients in
    iid_randomized mode.
    """
    if (min(sample_size, target_size, source_count) < 2
            or any(x != int(x) for x in [sample_size, target_size, source_count])
            or not 0 < valid_fraction < valid_probability < 1 or not 0 < alpha < .5
            or not np.isfinite(shift_constant) or shift_constant <= 0):
        raise ValueError('Valid sizes, mixture slack and positive shift constant required')
    if experiment not in ['fixed_arm','iid_randomized']:
        raise ValueError('Choose fixed_arm or iid_randomized experiment')
    factor = 4. if experiment=='fixed_arm' else 2.
    rate = min(1/np.sqrt(sample_size)/source_count**.25, 1/np.sqrt(target_size))
    shift = shift_constant*rate
    chi_square = mixture_chi_square(sample_size, shift, valid_probability, factor)
    log_joint_moment = source_count*np.log1p(chi_square)+target_size*np.log1p(factor*shift**2)
    total_variation_bound = 1. if log_joint_moment >= np.log(5.) else .5*np.sqrt(np.expm1(log_joint_moment))
    minimum_valid = int(np.ceil(valid_fraction*source_count))
    count_failure = float(binom.cdf(minimum_valid-1, source_count, valid_probability))
    bound = shift*max(0., 1-2*alpha-total_variation_bound-count_failure)
    return dict(source_size=sample_size, target_size=target_size, source_count=source_count,
        experiment=experiment, sample_size_units='treated patients' if factor==4 else 'all patients',
        shift_constant=shift_constant, shift=shift, maximum_source_bias=shift/(1-valid_probability),
        site_chi_square=chi_square, joint_log_second_moment=float(log_joint_moment),
        total_variation_bound=float(total_variation_bound), count_failure=count_failure,
        length_lower_bound=float(bound), scaled_lower_bound=float(bound/rate))


def draw_valid_flags(rng, source_count, valid_probability=.8, valid_fraction=.75):
    """Draw the count-conditioned lower-bound prior; its labels are audit-only."""
    if (source_count < 1 or source_count != int(source_count)
            or not 0 < valid_fraction < valid_probability < 1):
        raise ValueError('A positive integer count and strict valid-probability slack required')
    minimum_valid = int(np.ceil(valid_fraction*source_count))
    while True:
        flags = rng.uniform(size=source_count) < valid_probability
        if flags.sum() >= minimum_valid:
            return flags

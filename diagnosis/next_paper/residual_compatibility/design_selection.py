"""Outcome-independent site removal for conditional balanced inference.

Arm order is treated, control throughout. Original valid-source guarantees
are declared before design selection; removal never estimates site validity.
"""
import numpy as np

from continuous_balance import (UnsupportedDesign, linear_balance_design,
    evaluate_linear_balance, conditional_linear_bands)
from population_protection import projection_geometry, projection_error_budget


def select_source_designs(source_designs, target_mean, valid_minimum,
                          design_policy='all_required'):
    """Prepare both arms using only (feature matrix, treatment) site tuples.

    Each feature matrix already includes its intercept. A site is removed if
    either arm fails the existing rank, condition-number or leverage guard.
    Supported designs are cached even under all_required for paired audits.
    """
    if design_policy not in ['all_required', 'remove']:
        raise ValueError('design_policy must be all_required or remove')
    retained, rejected, prepared = [], [], []
    count = 0
    for index, (matrix, treatment) in enumerate(source_designs):
        count += 1
        matrix, treatment = np.asarray(matrix), np.asarray(treatment)
        if (matrix.ndim != 2 or treatment.shape != (len(matrix),)
                or not np.all(np.isin(treatment, [0, 1]))):
            raise ValueError('Aligned feature matrices and binary treatment required')
        arms = []
        try:
            for arm in [1, 0]:
                rows = np.flatnonzero(treatment == arm)
                arms.append(dict(rows=rows,
                    design=linear_balance_design(matrix[rows], target_mean)))
        except UnsupportedDesign as error:
            rejected.append(dict(index=index, reason=str(error)))
            continue
        retained.append(index)
        prepared.append(arms)
    valid = np.broadcast_to(valid_minimum, (2,)).astype(float)
    if (count < 1 or not np.all(np.isfinite(valid)) or np.any(valid < 1)
            or np.any(valid > count) or np.any(valid != np.floor(valid))):
        raise ValueError('Original integer valid-source counts in [1,K] required')
    adjusted = np.maximum(0, valid.astype(int)-len(rejected))
    eligible = bool(retained and np.all(adjusted > 0)
                    and (design_policy == 'remove' or not rejected))
    return dict(original_count=count, original_valid=valid.astype(int),
        retained_indices=np.array(retained, dtype=int), rejected=rejected,
        designs=prepared, retained_count=len(retained), adjusted_valid=adjusted,
        eligible=eligible, design_policy=design_policy)


def evaluate_retained_sources(selection, source_outcomes):
    """Evaluate cached designs after the entire selection has been fixed."""
    if len(source_outcomes) != selection['original_count']:
        raise ValueError('One outcome array per original site is required')
    if selection['retained_count'] == 0:
        return None
    site_summaries = []
    for index, arms in zip(selection['retained_indices'], selection['designs']):
        outcomes = np.asarray(source_outcomes[index])
        if outcomes.shape != (sum(len(arm['rows']) for arm in arms),):
            raise ValueError('Outcomes must align with the original site design')
        values = [evaluate_linear_balance(arm['design'], outcomes[arm['rows']])
                  for arm in arms]
        site_summaries.append({key: np.array([v[key] for v in values]) for key in values[0]})
    return {key: np.stack([value[key] for value in site_summaries])
            for key in site_summaries[0]}


def retained_projection_budgets(selection, source, model_mode='linear',
                                projection_envelopes=None):
    """Recompute unknown-invalid-subset certificates on the retained sites.

    Target and valid-source means remain linear in all three modes. Declared
    mode assumes |(I-H)m|_2 <= envelope_a sqrt(n_ja) for each invalid source.
    Bounded-invalid mode uses the universal envelope 1/2 for centered means.
    """
    if model_mode not in ['linear', 'declared', 'bounded_invalid']:
        raise ValueError('Unknown model mode')
    if model_mode == 'linear':
        if projection_envelopes is not None:
            raise ValueError('Projection envelopes apply only to declared mode')
        return None
    if model_mode == 'bounded_invalid':
        if projection_envelopes is not None:
            raise ValueError('Bounded-invalid mode does not use supplied envelopes')
        envelopes = np.array([.5, .5])
    else:
        if projection_envelopes is None:
            raise ValueError('Declared mode requires two projection envelopes')
        envelopes = np.broadcast_to(projection_envelopes, (2,)).astype(float)
        if not np.all(np.isfinite(envelopes)) or np.any(envelopes < 0):
            raise ValueError('Finite nonnegative projection envelopes required')
    budgets = {}
    for arm, name in enumerate(['mu1', 'mu0']):
        summaries = {key: value[:, arm] for key, value in source.items()}
        geometry = [projection_geometry(arms[arm]['design']) for arms in selection['designs']]
        for key in ['sample_size', 'projection_commutator_norm']:
            summaries[key] = np.array([item[key] for item in geometry])
        invalid_maximum = selection['retained_count']-selection['adjusted_valid'][arm]
        budgets[name] = projection_error_budget(summaries,
            envelopes[arm]*np.sqrt(summaries['sample_size']), invalid_maximum)
    budgets['contrast'] = {key: budgets['mu1'][key]+budgets['mu0'][key]
                           for key in budgets['mu1']}
    return budgets


def conditional_selected_bands(selection, source, target, model_mode='linear',
                                projection_envelopes=None, source_failure=.01,
                                target_failure=.01, dispersion_failure=.015):
    """Recompute all normalizations and dispersion geometry after removal."""
    if model_mode not in ['linear', 'declared', 'bounded_invalid']:
        raise ValueError('Unknown model mode')
    if not selection['eligible']:
        if not 0 < target_failure < 1:
            raise ValueError('Target failure allowance must be in (0,1)')
        point = float(target['estimate'][0]-target['estimate'][1])
        noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/target_failure))
        band = point+noise*np.array([-1., 1.])
        return dict(conditional_interval=band, target_band=band,
            source_band=np.array([-np.inf, np.inf]), source_noise=0.,
            dispersion_radius=0., target_projection=point, target_point=point,
            target_noise=noise, fallback=False, design_fallback=True)
    if source is None or source['estimate'].shape != (selection['retained_count'], 2):
        raise ValueError('Source summaries must match the retained site set')
    budgets = retained_projection_budgets(selection, source, model_mode, projection_envelopes)
    result = conditional_linear_bands(source, target, selection['adjusted_valid'],
        source_failure=source_failure, target_failure=target_failure,
        dispersion_failure=dispersion_failure, dispersion_budgets=budgets)
    result['design_fallback'] = False
    return result

"""Population inference with a common composition event across design branches."""
import numpy as np

from design_selection import conditional_selected_bands
from population_protection import composition_radius, expand_conditional_interval


def inference_ledger(eligible, budget_policy='shared_composition'):
    """Prespecified 5% allocations; conditional errors can depend on design.

    shared_composition gives .005 to the coefficient event and .0225 to a
    common composition event. The remaining .0225 is conditional on design.
    Its target fallback is exactly the standalone protected target procedure.
    """
    if budget_policy not in ['legacy', 'conditional_reallocation', 'shared_composition']:
        raise ValueError('Unknown error-budget policy')
    composition = .0225 if budget_policy == 'shared_composition' else .01
    conditional = .0225 if budget_policy == 'shared_composition' else .035
    scale = conditional/.035
    source, target, dispersion = np.array([.01, .01, .015])*scale
    if not eligible and budget_policy != 'legacy':
        target = conditional
        source, dispersion = 0., 0.
    return dict(coefficient=.005, composition=composition, source=float(source),
        target=float(target), dispersion=float(dispersion), conditional=conditional)


def design_summary(selection):
    """Set observed dispersion to zero while preserving only design factors."""
    fields = ['weight_square_sum', 'diagonal_square_sum', 'correction_frobenius',
              'max_diagonal', 'balance_error', 'negative_weight_fraction']
    shape = (selection['retained_count'], 2)
    summary = {name: np.array([[arm['design'][name] for arm in arms]
                              for arms in selection['designs']]) for name in fields}
    summary.update(estimate=np.zeros(shape), variance=np.zeros(shape))
    return summary


def apply_design_precision_gate(selection, target_weight_square_sum,
                                model_mode='linear', projection_envelopes=None,
                                budget_policy='shared_composition'):
    """Compare zero-observed-dispersion reference radius with target noise.

    This is a design-only heuristic for potential precision, not a lower bound
    on the realized interval width. It cannot guarantee pointwise noninferiority.
    """
    if not selection['eligible']:
        return dict(selection, precision_gate_rejected=False)
    target_weights = np.asarray(target_weight_square_sum, dtype=float)
    if target_weights.shape != (2,) or np.any(target_weights <= 0) or not np.all(np.isfinite(target_weights)):
        raise ValueError('Two finite positive target weight norms required')
    ledger = inference_ledger(True, budget_policy)
    target = dict(estimate=np.zeros(2), weight_square_sum=target_weights)
    reference = conditional_selected_bands(selection, design_summary(selection), target,
        model_mode=model_mode, projection_envelopes=projection_envelopes,
        source_failure=ledger['source'], target_failure=ledger['target'], dispersion_failure=ledger['dispersion'])
    source_radius = reference['source_noise']+reference['dispersion_radius']
    target_radius = np.sqrt(target_weights.sum()/2*np.log(2/ledger['conditional']))
    rejected = bool(source_radius >= target_radius)
    return dict(selection, eligible=not rejected, precision_gate_rejected=rejected,
        design_source_radius=float(source_radius), design_target_radius=float(target_radius))


def population_selected_interval(selection, source, target, certificate, target_size,
                                  model_mode='linear', projection_envelopes=None,
                                  budget_policy='shared_composition'):
    """Use the same unconditional composition budget on both design branches."""
    ledger = inference_ledger(selection['eligible'], budget_policy)
    if not np.isclose(certificate['failure'], ledger['coefficient'], rtol=0, atol=1e-14):
        raise ValueError('Coefficient certificate must use the declared ledger')
    result = conditional_selected_bands(selection, source, target, model_mode=model_mode,
        projection_envelopes=projection_envelopes, source_failure=ledger['source'],
        target_failure=ledger['target'], dispersion_failure=ledger['dispersion'])
    radius = composition_radius(certificate, target_size, ledger['composition'])
    interval = expand_conditional_interval(result['conditional_interval'], radius)
    return dict(result, interval=interval, composition_radius=radius, ledger=ledger,
        point=float(np.clip(result['target_projection'], interval[0], interval[1])))

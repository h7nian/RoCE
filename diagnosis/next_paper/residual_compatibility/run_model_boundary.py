"""Exact contamination-overlap audit and paired finite-support model comparison.

The lower-bound audit is conditional on distinct fixed patient covariates.
The saturated-basis comparison assumes only three covariate support points.
These are different model classes; neither establishes all-valid optimality.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import binom

from continuous_balance import linear_balance_design, conditional_linear_bands
from population_protection import projection_geometry, projection_error_budget
from run_projection_stress import fixed_geometry, summaries


def contamination_overlap(arm_size, target_arm_size, source_count, epsilon=.2,
                           valid_fraction=.75, shift_constant=.1, alpha=.05):
    """Finite binomial calculation for the two-prior interval-length bound."""
    if (min(arm_size, target_arm_size, source_count) < 1 or not 0 < epsilon < 1
            or not 0 < valid_fraction < 1-epsilon or not 0 < alpha < .5
            or not 0 < shift_constant < .5):
        raise ValueError('Positive sizes and a strict valid-fraction slack required')
    delta = shift_constant/np.sqrt(max(arm_size, target_arm_size))
    support = np.arange(arm_size+1)
    p0, p1 = [binom.pmf(support, arm_size, probability) for probability in [.5, .5+delta]]
    # Remove floating-point summation error in scipy's log-gamma PMF evaluation.
    p0, p1 = p0/p0.sum(), p1/p1.sum()
    total_variation = float(np.abs(p1-p0).sum()/2)
    slack = 1-(1-epsilon)*(1+total_variation)
    if slack < -1e-12:
        raise ValueError('The contamination neighborhoods do not overlap')
    common = (1-epsilon)*np.maximum(p0, p1)+max(0., slack)*p0
    q0, q1 = (common-(1-epsilon)*p0)/epsilon, (common-(1-epsilon)*p1)/epsilon
    target_support = np.arange(target_arm_size+1)
    target_tv = float(np.abs(binom.pmf(target_support, target_arm_size, .5+delta)
                           -binom.pmf(target_support, target_arm_size, .5)).sum()/2)
    minimum_valid = int(np.ceil(valid_fraction*source_count))
    count_failure = float(binom.cdf(minimum_valid-1, source_count, 1-epsilon))
    length_bound = delta*max(0., 1-2*alpha-target_tv-2*count_failure)
    return dict(source_arm_size=arm_size, target_arm_size=target_arm_size,
        source_count=source_count, delta=delta, source_tv=total_variation,
        target_tv=target_tv, count_failure=count_failure,
        interval_length_lower_bound=length_bound,
        scaled_lower_bound=length_bound*np.sqrt(max(arm_size, target_arm_size)),
        common_mass=float(common.sum()), q0_mass=float(q0.sum()), q1_mass=float(q1.sum()),
        minimum_contamination_mass=float(min(q0.min(), q1.min())),
        mixture_error=float(np.max(np.abs((1-epsilon)*p0+epsilon*q0-((1-epsilon)*p1+epsilon*q1)))))


def saturated_geometry(arm_size=500):
    values, counts, _, _ = fixed_geometry(arm_size)
    x = np.repeat(values, counts)
    matrix = np.column_stack((np.ones(arm_size), x, x*x))
    design = linear_balance_design(matrix, matrix.mean(axis=0))
    diagonal = np.array([design['diagonal'][np.flatnonzero(x == value)[0]] for value in values])
    constants = {name: design[name] for name in ['weight_square_sum', 'diagonal_square_sum',
        'correction_frobenius', 'max_diagonal', 'balance_error', 'negative_weight_fraction']}
    constants.update(projection_geometry(design))
    return values, counts, diagonal, constants


def saturated_summaries(successes, counts, diagonal, constants):
    mean = successes.sum(axis=-1)/counts.sum()
    variance = np.sum(diagonal*successes*(1-successes/counts), axis=-1)
    return dict(estimate=mean, variance=variance,
        **{name: np.broadcast_to(value, mean.shape) for name, value in constants.items()})


def run_basis_cell(count, eta, repeats, seed):
    values, counts, diagonal, constants = fixed_geometry()
    _, _, saturated_diagonal, saturated_constants = saturated_geometry()
    template = {name: np.full(count, value) for name, value in constants.items()}
    template['estimate'] = np.zeros(count)
    budget = projection_error_budget(template, .5*np.sqrt(counts.sum()), count//4)
    budgets = dict(mu1=budget, mu0=budget,
        contrast={key: 2*value for key, value in budget.items()})
    probabilities = np.full((count, 3), .5)
    probabilities[:count//4] += 2*np.sqrt(.5/1000)*count**(-.25)+eta*(values**2-.5)
    records = []
    for repeat in range(repeats):
        rng = np.random.RandomState(seed+repeat)
        target_draws = [rng.binomial(counts, .5), rng.binomial(counts, .4)]
        source_draws = [rng.binomial(counts, probabilities), rng.binomial(counts, np.full((count, 3), .4))]
        variants = []
        for saturated in [False, True]:
            def summarize(draw):
                if saturated:
                    return saturated_summaries(draw, counts, saturated_diagonal, saturated_constants)
                return summaries(draw, values, counts, diagonal, constants)
            target_arms = [summarize(draw) for draw in target_draws]
            source_arms = [summarize(draw) for draw in source_draws]
            target = {key: np.array([arm[key] for arm in target_arms]) for key in target_arms[0]}
            source = {key: np.stack([arm[key] for arm in source_arms], axis=-1) for key in source_arms[0]}
            result = conditional_linear_bands(source, target, 3*count//4,
                dispersion_budgets=None if saturated else budgets)
            variants.append({key: result[key] for key in ['conditional_interval', 'source_band', 'dispersion_radius']})
        records.append(variants)
    arrays = {key: np.array([[v[key] for v in record] for record in records])
              for key in ['conditional_interval', 'source_band', 'dispersion_radius']}
    metrics = []
    for j, label in enumerate(['linear_basis_bounded_invalid', 'saturated_basis']):
        interval, source = arrays['conditional_interval'][:, j], arrays['source_band'][:, j]
        metrics.append(dict(source_count=count, nonlinearity=eta, method=label, repeats=repeats,
            conditional_coverage=float(np.mean((interval[:, 0] <= .1) & (.1 <= interval[:, 1]))),
            source_coverage=float(np.mean((source[:, 0] <= .1) & (.1 <= source[:, 1]))),
            conditional_length=float(np.diff(interval, axis=1).mean()), source_length=float(np.diff(source, axis=1).mean())))
    return arrays, metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--repeats', type=int, default=200)
    args = parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats < 2:
        raise ValueError('Use FACE-HD scratch and >=2 repetitions')
    output = args.root/'model_boundary'
    output.mkdir()
    counts = [128, 2048, 32768, 262144]
    (output/'configuration.json').write_text(json.dumps(dict(repeats=args.repeats, source_counts=counts,
        amplitudes=[0., .4], seed_offset=32190000, arm_size=500,
        comparison='Same patients; basis fixed in advance; three-point support, not unrestricted continuous means'), indent=2)+'\n')
    lower_bounds = [contamination_overlap(n, n, k) for n in [125, 500, 2000] for k in counts]
    (output/'contamination_overlap.json').write_text(json.dumps(lower_bounds, indent=2)+'\n')
    rows = []
    for index, count in enumerate(counts):
        for eta in [0., .4]:
            arrays, metrics = run_basis_cell(count, eta, args.repeats, 32190000+10000*index)
            stem = output/('k{}_eta{}'.format(count, int(100*eta)))
            np.savez_compressed(str(stem)+'.npz', **arrays)
            Path(str(stem)+'.json').write_text(json.dumps(dict(metrics=metrics), indent=2)+'\n')
            rows.extend(metrics)
            print('Boundary comparison', count, eta, flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (output/'COMPLETE').write_text('Conditional coverage only; no population enlargement conceals source failures.\n')


if __name__ == '__main__':
    main()

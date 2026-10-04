"""Fresh-design continuous-covariate checks for the linear-balance benchmark."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from scipy.stats import norm

from continuous_balance import (UnsupportedDesign, linear_balance_design,
    evaluate_linear_balance, continuous_balance_interval, leaveout_dispersion_upper)


SUMMARY_FIELDS = ('estimate', 'variance', 'weight_square_sum', 'diagonal_square_sum',
                  'correction_frobenius', 'max_diagonal', 'balance_error', 'negative_weight_fraction')


def stack_summaries(arms):
    return {name: np.array([arm[name] for arm in arms]) for name in SUMMARY_FIELDS}


def draw_site(rng, setting, target=False, index=0):
    size, dimension = setting['n_per_site'], setting['dimension']
    features = rng.uniform(-1., 1., (size, dimension))
    if not target:
        shifts = setting.get('source_shift', .4) * np.array([1., -.6, .4, -.3])
        for coordinate, shift in enumerate(shifts):
            if shift != 0:
                uniform = (features[:, coordinate] + 1) / 2
                features[:, coordinate] = (-1 + np.sqrt((1-shift)**2 + 4*shift*uniform)) / shift
    propensity = .5 + .12*features[:, 0] - .08*features[:, 1]
    if not target:
        propensity += .04*np.cos(index)
    treatment = rng.binomial(1, propensity)
    base = .35 + features[:, :4].dot(np.array([.06, -.035, .025, .02]))
    probabilities = np.column_stack((base + .1 + setting['heterogeneity']*features[:, 0], base))
    if setting['pattern'] == 'nonlinear':
        probabilities += .15*features[:, 0, None]*features[:, 1, None]
    if not target:
        count = setting['source_count']
        quarter = count // 4
        weak = 2*np.sqrt(.5/size)*count**(-.25)
        pattern = setting['pattern']
        if pattern == 'weak' and index < quarter:
            probabilities[:, 0] += weak
        elif pattern == 'strong' and index < quarter:
            probabilities[:, 0] += .15
        elif pattern == 'cancel' and index < quarter:
            probabilities += .15
        elif pattern == 'disjoint':
            if index < quarter:
                probabilities[:, 0] += weak
            elif index < 2*quarter:
                probabilities[:, 1] -= weak
    if np.any((probabilities < 0) | (probabilities > 1)):
        raise ArithmeticError('DGP produced invalid outcome probabilities')
    observed = rng.binomial(1, probabilities[np.arange(size), 1-treatment])
    return features, treatment, observed, probabilities


def fit_site(site, target_mean):
    features, treatment, outcomes, probabilities = site
    matrix = np.column_stack((np.ones(len(features)), features))
    summaries, coefficients, conditional_means, conditional_variances = [], [], [], []
    leverage = 0.
    for arm in [1, 0]:
        rows = treatment == arm
        design = linear_balance_design(matrix[rows], target_mean)
        summaries.append(evaluate_linear_balance(design, outcomes[rows]))
        coefficients.append(np.linalg.solve(design['triangular'], design['basis'].T.dot(outcomes[rows])))
        means = probabilities[rows, 1-arm]
        conditional_means.append(design['weights'].dot(means))
        conditional_variances.append(np.sum(design['weights']**2*means*(1-means)))
        leverage = max(leverage, design['leverage'].max())
    # The last two items are oracle diagnostics and are never passed to inference.
    return stack_summaries(summaries), np.array(coefficients), np.array(conditional_means), np.array(conditional_variances), leverage


def fallback_interval(target, target_size, cate_range):
    point = target['estimate'][0] - target['estimate'][1]
    noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(200.))
    composition = cate_range*np.sqrt(np.log(200.)/(2*target_size))
    interval = np.clip(point + (noise+composition)*np.array([-1., 1.]), -1, 1)
    return interval, point


def run_repeat(setting, repeat):
    rng = np.random.RandomState(setting['seed_offset'] + 10000*setting['cell_id'] + repeat)
    target_site = draw_site(rng, setting, target=True)
    target_features = target_site[0]
    target_mean = np.r_[1., target_features.mean(axis=0)]
    conditional_truth = np.mean(target_site[3][:, 0] - target_site[3][:, 1])
    try:
        target, coefficients, _, target_variance, max_leverage = fit_site(target_site, target_mean)
    except UnsupportedDesign:
        # Design-only fallback; keep this repetition, including in every summary.
        methods = {name: dict(interval=np.array([-1., 1.]), point=0., midpoint=0.,
                   design_fallback=True, conditional_coverage=True)
                   for name in ['target_universal', 'arms_universal', 'adaptive_universal',
                                'target_declared_half', 'arms_declared_half',
                                'adaptive_declared_half', 'target_normal']}
        return dict(methods=methods, diagnostics=dict(conditional_truth=conditional_truth,
            target_failed=True, source_failed=True, naive_conditional_coverage=True))
    fitted_contrast = np.column_stack((np.ones(len(target_features)), target_features)).dot(coefficients[0]-coefficients[1])
    target_variance_estimate = max(0., target['variance'].sum() + fitted_contrast.var(ddof=1)/setting['n_per_site'])
    target_point = target['estimate'][0]-target['estimate'][1]
    target_normal = np.clip(target_point + norm.isf(.025)*np.sqrt(target_variance_estimate)*np.array([-1., 1.]), -1, 1)
    values, means, variances = [], [], []
    source_failed = False
    for index in range(setting['source_count']):
        site = draw_site(rng, setting, index=index)
        try:
            summary, _, expectation, variance, leverage = fit_site(site, target_mean)
        except UnsupportedDesign:
            source_failed = True
            continue
        values.append(summary); means.append(expectation); variances.append(variance)
        max_leverage = max(max_leverage, leverage)
    methods = {}
    diagnostics = dict(conditional_truth=conditional_truth, max_leverage=max_leverage,
                       target_variance_estimate=target_variance_estimate,
                       target_variance_oracle=float(target_variance.sum()), source_failed=source_failed, target_failed=False)
    for bound in [2., .5]:
        label = 'universal' if bound == 2. else 'declared_half'
        baseline_noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(80.))
        baseline_composition = bound*np.sqrt(np.log(80.)/(2*setting['n_per_site']))
        baseline = np.clip(target_point + (baseline_noise+baseline_composition)*np.array([-1., 1.]), -1, 1)
        methods['target_'+label] = dict(interval=baseline, point=target_point)
        for mode in ['arms', 'adaptive']:
            name = mode+'_'+label
            if source_failed:
                interval, point = fallback_interval(target, setting['n_per_site'], bound)
                methods[name] = dict(interval=interval, point=point, midpoint=interval.mean(),
                                     design_fallback=True, conditional_coverage=True)
                continue
            source = {key: np.stack([value[key] for value in values]) for key in SUMMARY_FIELDS}
            result = continuous_balance_interval(source, target, 3*setting['source_count']//4,
                setting['n_per_site'], dispersion_mode=mode, cate_range=bound)
            methods[name] = dict(interval=result['interval'], point=result['target_projection'],
                midpoint=result['midpoint'], design_fallback=False,
                source_noise=result['source_noise'], dispersion_radius=result['dispersion_radius'],
                conditional_coverage=bool(result['source_band'][0] <= conditional_truth <= result['source_band'][1]))
    methods['target_normal'] = dict(interval=target_normal, point=target_point)
    if source_failed:
        diagnostics['naive_conditional_coverage'] = True  # Explicit full-range diagnostic fallback.
    if not source_failed:
        expected = np.array(means)
        source = {key: np.stack([value[key] for value in values]) for key in SUMMARY_FIELDS}
        source_mean = np.mean(source['estimate'][:, 0]-source['estimate'][:, 1])
        # Conditional-only naive source interval: label distinguishes its estimand.
        source_variance = max(0., source['variance'].sum()/setting['source_count']**2)
        naive = source_mean + norm.isf(.025)*np.sqrt(source_variance)*np.array([-1., 1.])
        diagnostics['naive_conditional_coverage'] = bool(naive[0] <= conditional_truth <= naive[1])
        diagnostics['source_variance_estimate'] = source_variance
        diagnostics['source_variance_oracle'] = float(np.sum(variances)/setting['source_count']**2)
        diagnostics['maximum_balance_error'] = float(source['balance_error'].max())
        diagnostics['negative_weight_fraction'] = float(source['negative_weight_fraction'].mean())
        diagnostics['dispersion_failures'] = 0
        for arm in range(2):
            upper = leaveout_dispersion_upper({key: value[:, arm] for key, value in source.items()}, .005)
            diagnostics['dispersion_failures'] += int(np.var(expected[:, arm]) > upper)
        guaranteed = np.arange(setting['source_count']) >= setting['source_count']//2
        target_arm_truth = target_site[3].mean(axis=0)
        diagnostics['maximum_valid_drift'] = float(np.max(np.abs(expected[guaranteed]-target_arm_truth)))
    return dict(methods=methods, diagnostics=diagnostics)


def run_cell(task):
    setting, repeats, output = task
    started = time.time()
    records = [run_repeat(setting, repeat) for repeat in range(repeats)]
    prefix = Path(output)/('cell_{:03d}'.format(setting['cell_id']))
    names = list(records[0]['methods'])
    arrays = dict(intervals=np.array([[r['methods'][name]['interval'] for name in names] for r in records]),
                  points=np.array([[r['methods'][name]['point'] for name in names] for r in records]),
                  midpoints=np.array([[r['methods'][name].get('midpoint', r['methods'][name]['point']) for name in names] for r in records]))
    np.savez_compressed(str(prefix)+'.npz', **arrays)
    metrics = []
    for index, name in enumerate(names):
        intervals = arrays['intervals'][:, index]
        errors = arrays['points'][:, index]-.1
        coverage = np.mean((intervals[:, 0] <= .1) & (.1 <= intervals[:, 1]))
        metrics.append(dict(setting, method=name, repeats=repeats, coverage=float(coverage),
            coverage_mcse=float(np.sqrt(coverage*(1-coverage)/repeats)),
            length=float(np.mean(intervals[:, 1]-intervals[:, 0])), bias=float(errors.mean()),
            rmse=float(np.sqrt(np.mean(errors**2))),
            midpoint_rmse=float(np.sqrt(np.mean((arrays['midpoints'][:, index]-.1)**2)))))
    Path(str(prefix)+'.json').write_text(json.dumps(dict(setting=setting, methods=names, metrics=metrics,
        diagnostics=[r['diagnostics'] for r in records], runtime_seconds=time.time()-started), indent=2)+'\n')
    return metrics


def plan(profile):
    if profile == 'main':
        cells = [(p, k, pattern, .15, .4) for p in [4, 10, 20]
                 for k in [16, 64, 256] for pattern in ['all_valid', 'weak']]
    elif profile == 'boundary':
        cells = [(p, 64, 'weak', .15, .4) for p in [50, 100, 200]]
        cells += [(20, 64, pattern, .15, .4) for pattern in ['strong', 'cancel', 'disjoint', 'nonlinear']]
        cells += [(20, 64, 'weak', .15, shift) for shift in [0., .8]]
    elif profile == 'pilot':
        cells = [(4, 16, 'weak', .15, .4), (20, 64, 'weak', .15, .4), (200, 64, 'weak', .15, .4)]
    elif profile == 'stress':
        cells = [(10, k, pattern, 0., .4) for k in [128, 512, 2048] for pattern in ['all_valid', 'weak']]
    else:
        raise ValueError('Unknown experiment profile')
    return [dict(cell_id=i+1, dimension=p, source_count=k, pattern=pattern,
                 heterogeneity=h, source_shift=shift, n_per_site=1000,
                 seed_offset=dict(main=19190000, boundary=20190000, pilot=21190000, stress=22190000)[profile])
            for i, (p, k, pattern, h, shift) in enumerate(cells)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--profile', choices=['pilot', 'main', 'boundary', 'stress'], default='main')
    parser.add_argument('--repeats', type=int, default=200)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats < 2 or not 1 <= args.workers <= 4:
        raise ValueError('Use FACE-HD scratch, >=2 repeats and 1--4 workers')
    output = args.root/args.profile
    output.mkdir()
    settings = plan(args.profile)
    config = dict(settings=settings, repeats=args.repeats,
        scope='All-site linear conditional-mean benchmark, not sparse high-dimensional or nonlinear inference',
        population_truth=.1, declared_cate_ranges=[2., .5],
        code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in Path(__file__).parent.glob('*continuous_balance*.py')})
    (output/'configuration.json').write_text(json.dumps(config, indent=2)+'\n')
    rows = []
    with multiprocessing.Pool(args.workers) as pool:
        for metrics in pool.imap_unordered(run_cell, [(s, args.repeats, str(output)) for s in settings]):
            rows.extend(metrics)
            print('Completed', args.profile, metrics[0]['cell_id'], flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(sorted(rows, key=lambda r: (r['cell_id'], r['method'])))
    (output/'COMPLETE').write_text('All repetitions retained.\n')


if __name__ == '__main__':
    main()

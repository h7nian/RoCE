"""Independent confirmation and evaluated design-removal experiments."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from scipy.special import expit
from scipy.stats import norm

from continuous_balance import UnsupportedDesign
from design_selection import (select_source_designs, evaluate_retained_sources,
    conditional_selected_bands)
from population_protection import target_coefficient_certificate, composition_radius, expand_conditional_interval
from run_target_protection import fit_packet, run_cell as run_confirmation_cell


def draw_design_site(rng, setting, target=False, index=0):
    size = 1000 if target else setting['source_size']
    features = rng.uniform(-1, 1, (size, setting['dimension']))
    if not target and index >= setting['source_count']*(1-setting.get('singular_fraction', 0.)):
        # A prespecified source subpopulation with a rank-deficient feature map.
        features[:, -1] = features[:, -2]
    strength = 0. if target else setting['overlap_strength']
    treatment = rng.binomial(1, expit(strength*(features[:, 0]+features[:, 1])))
    base = .35+features[:, :4].dot([.06, -.035, .025, .02])
    probabilities = np.column_stack((base+.1+setting['heterogeneity']*features[:, 0], base))
    if not target and setting['pattern'] == 'weak' and index < setting['source_count']//4:
        probabilities[:, 0] += 2*np.sqrt(.5/size)*setting['source_count']**(-.25)
    if np.any((probabilities < 0) | (probabilities > 1)):
        raise ArithmeticError('DGP outcome probabilities must stay in [0,1]')
    outcomes = rng.binomial(1, probabilities[np.arange(size), 1-treatment])
    return features, treatment, outcomes, probabilities


def run_removal_repeat(setting, repeat, design_policies, model_modes):
    rng = np.random.RandomState(31190000+10000*setting['cell_id']+repeat)
    target_site = draw_design_site(rng, setting, target=True)
    features = target_site[0]
    target_mean = np.r_[1., features.mean(axis=0)]
    truth = float(np.mean(target_site[3][:, 0]-target_site[3][:, 1]))
    labels = [policy+'_'+mode for policy in design_policies for mode in model_modes]
    labels += ['target_matched', 'target_protected', 'target_normal']
    try:
        target, _, _ = fit_packet(target_site, target_mean)
        certificate = target_coefficient_certificate(features, target_site[1], target_site[2])
    except UnsupportedDesign:
        methods = {name: dict(interval=[-1., 1.], conditional=[-1., 1.], point=0.,
            source_noise=0., dispersion=0., composition=1., design_fallback=True) for name in labels}
        return dict(methods=methods, diagnostics=dict(conditional_truth=truth, target_failed=True,
            rejected_count=setting['source_count'], retained_count=0, adjusted_valid=[0, 0],
            eligible=False, retained_indices=[], max_weight_square_sum=None))
    sites = [draw_design_site(rng, setting, index=j) for j in range(setting['source_count'])]
    # This function accepts no outcomes or DGP probabilities.
    selection = select_source_designs(
        ((np.column_stack((np.ones(len(site[0])), site[0])), site[1]) for site in sites),
        target_mean, 3*setting['source_count']//4, design_policy='remove')
    source = evaluate_retained_sources(selection, [site[2] for site in sites])
    methods = {}
    radius = composition_radius(certificate, len(features), .01)
    for policy in design_policies:
        chosen = dict(selection, design_policy=policy,
            eligible=selection['eligible'] and (policy == 'remove' or not selection['rejected']))
        for mode in model_modes:
            result = conditional_selected_bands(chosen, source, target, model_mode=mode)
            conditional = result['conditional_interval']
            interval = expand_conditional_interval(conditional, radius)
            methods[policy+'_'+mode] = dict(interval=interval,
                conditional=conditional, point=float(np.clip(result['target_projection'], interval[0], interval[1])),
                source_noise=float(result['source_noise']), dispersion=float(result['dispersion_radius']),
                composition=radius, design_fallback=result['design_fallback'])
    point = float(target['estimate'][0]-target['estimate'][1])
    for name, failure in [('target_matched', .01), ('target_protected', .0225)]:
        noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/failure))
        conditional = point+noise*np.array([-1., 1.])
        composition = composition_radius(certificate, len(features), failure)
        methods[name] = dict(interval=expand_conditional_interval(conditional, composition),
            conditional=conditional, point=point, source_noise=0., dispersion=0.,
            composition=composition, design_fallback=False)
    variance = max(0., target['variance'].sum()+np.var(features.dot(certificate['coefficient']), ddof=1)/len(features))
    methods['target_normal'] = dict(interval=np.clip(point+norm.isf(.025)*np.sqrt(variance)*np.array([-1., 1.]), -1, 1),
        conditional=point+norm.isf(.025)*np.sqrt(max(0., target['variance'].sum()))*np.array([-1., 1.]),
        point=point, source_noise=0., dispersion=0., composition=0., design_fallback=False)
    audit = dict(conditional_truth=truth, target_failed=False,
        rejected_count=len(selection['rejected']), retained_count=selection['retained_count'],
        adjusted_valid=selection['adjusted_valid'].tolist(), eligible=selection['eligible'],
        retained_indices=selection['retained_indices'].tolist(),
        max_weight_square_sum=float(source['weight_square_sum'].max()) if source is not None else None)
    return dict(methods=methods, diagnostics=audit)


def run_removal_cell(task):
    setting, repeats, output, policies, modes = task
    started = time.time()
    records = [run_removal_repeat(setting, repeat, policies, modes) for repeat in range(repeats)]
    methods = list(records[0]['methods'])
    arrays = {name: np.array([[r['methods'][m][name] for m in methods] for r in records])
              for name in ['interval', 'conditional', 'point', 'source_noise', 'dispersion', 'composition', 'design_fallback']}
    arrays['conditional_truth'] = np.array([r['diagnostics']['conditional_truth'] for r in records])
    metrics = []
    for j, method in enumerate(methods):
        interval, conditional = arrays['interval'][:, j], arrays['conditional'][:, j]
        metrics.append(dict(setting, method=method, repeats=repeats,
            coverage=float(np.mean((interval[:, 0] <= .1) & (.1 <= interval[:, 1]))),
            conditional_coverage=float(np.mean((conditional[:, 0] <= arrays['conditional_truth'])
                & (arrays['conditional_truth'] <= conditional[:, 1]))),
            length=float(np.diff(interval, axis=1).mean()), conditional_length=float(np.diff(conditional, axis=1).mean()),
            rmse=float(np.sqrt(np.mean((arrays['point'][:, j]-.1)**2))),
            design_fallback=float(arrays['design_fallback'][:, j].mean())))
    stem = Path(output)/('cell_%03d'%setting['cell_id'])
    np.savez_compressed(str(stem)+'.npz', **arrays)
    Path(str(stem)+'.json').write_text(json.dumps(dict(setting=setting, methods=methods, metrics=metrics,
        diagnostics=[r['diagnostics'] for r in records], runtime_seconds=time.time()-started), indent=2)+'\n')
    return metrics


def plan(profile):
    if profile == 'confirmation':
        cells = [dict(dimension=4, source_count=2048, heterogeneity=0., pattern=pattern)
                 for pattern in ['weak', 'all_valid']]
        cells += [dict(dimension=10, source_count=2048, heterogeneity=0., pattern='weak'),
                  dict(dimension=4, source_count=256, heterogeneity=.15, pattern='weak')]
        return [dict(s, cell_id=i+1, n_per_site=1000, source_shift=.4, seed_offset=30190000)
                for i, s in enumerate(cells)]
    cells = [dict(source_size=120, dimension=45, source_count=k, overlap_strength=2.5,
                  heterogeneity=.15, pattern=pattern, singular_fraction=0.)
             for k in [4, 16, 64] for pattern in ['all_valid', 'weak']]
    cells += [dict(source_size=200, dimension=45, source_count=64, overlap_strength=2.5,
                   heterogeneity=.15, pattern=pattern, singular_fraction=0.) for pattern in ['all_valid', 'weak']]
    cells += [dict(source_size=1000, dimension=4, source_count=256, overlap_strength=0.,
                   heterogeneity=.15, pattern='weak', singular_fraction=f) for f in [0., .125]]
    return [dict(s, cell_id=i+1) for i, s in enumerate(cells)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--profile', choices=['confirmation', 'removal'], required=True)
    parser.add_argument('--repeats', type=int, default=200)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--run-name')
    parser.add_argument('--design-policy', choices=['all_required', 'remove', 'paired'], default='paired')
    parser.add_argument('--model-mode', choices=['linear', 'bounded_invalid', 'paired'], default='paired')
    args = parser.parse_args()
    if (Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents
            or args.repeats < 2 or not 1 <= args.workers <= 4):
        raise ValueError('Use FACE-HD scratch, >=2 repeats and 1--4 workers')
    name = args.run_name or args.profile
    if Path(name).name != name or name in ['.', '..']:
        raise ValueError('Use a new single output-directory name')
    output = args.root/name
    output.mkdir()
    settings = plan(args.profile)
    policies = ['all_required', 'remove'] if args.design_policy == 'paired' else [args.design_policy]
    modes = ['linear', 'bounded_invalid'] if args.model_mode == 'paired' else [args.model_mode]
    (output/'configuration.json').write_text(json.dumps(dict(settings=settings, repeats=args.repeats,
        design_policies=policies, model_modes=modes,
        code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}), indent=2)+'\n')
    tasks = [(s, args.repeats, str(output)) for s in settings]
    if args.profile == 'removal':
        tasks = [task+(policies, modes) for task in tasks]
    worker = run_confirmation_cell if args.profile == 'confirmation' else run_removal_cell
    rows = []
    with multiprocessing.Pool(args.workers) as pool:
        for result in pool.imap_unordered(worker, tasks):
            rows.extend(result)
            print('Completed', args.profile, result[0]['cell_id'], flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader(); writer.writerows(rows)
    (output/'COMPLETE').write_text('All repetitions, conditional intervals and design fallbacks retained.\n')


if __name__ == '__main__':
    main()

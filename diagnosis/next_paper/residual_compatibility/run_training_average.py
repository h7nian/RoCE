"""Paired fixed-valid-subset and conditional-outcome remainder experiments."""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
import numpy as np
from fitted_strata import fitted_stratum_data
from calibration_remainders import calibration_remainder_budgets
from conditional_remainders import conditional_remainder_budgets, full_target_conditional_interval, reuse_source_design_for_evaluation
from training_average import binomial_ratio_mgf, training_average_budgets
from repairs import infer_repair
from run_repairs import save_case


def make_plan():
    rows = []
    for count in [8, 32, 128, 512]:
        for pattern in ['all_valid', 'weak_boundary']:
            for floor in [0., .3]:
                rows.append(dict(cell_id=len(rows)+1, panel='fitted_training_average',
                    source_count=count, pattern=pattern, shared_scale=floor, n_per_site=1000,
                    source_training_size=500, source_evaluation_size=250, target_evaluation_size=500,
                    pilot_size=250, valid_fraction=.75, seed_offset=10130000,
                    certified_reference='full_target_conditional'))
    return rows


def run_cell(task):
    setting, repeats, output, envelopes = task
    started = time.time()
    data = fitted_stratum_data(setting, repeats)
    if setting.get("reuse_source_design", False):
        data = reuse_source_design_for_evaluation(data)
    valid_count = int(math.ceil(.75 * setting['source_count']))
    reference = full_target_conditional_interval(data)
    data['certified_reference_interval'] = reference
    data['certified_reference_width'] = np.diff(reference, axis=1)[:, 0]
    results, diagnostics = {}, {}
    configurations = [('previous_corner', calibration_remainder_budgets(data, valid_count))]
    configurations.append(('training_corner_p004', training_average_budgets(data, valid_count, envelopes['004'])))
    configurations.append(('conditional_pilot', conditional_remainder_budgets(data, valid_count)))
    for key in ['002', '004', '006']:
        configurations.append(('conditional_training_p'+key, conditional_remainder_budgets(data, valid_count,
            source_calibration='training_mgf', envelope=envelopes[key])))
    actual_mean = np.empty((repeats, 2))
    for arm in [0, 1]:
        indices = np.flatnonzero(data['structural_shift'][:, arm] == 0)[:valid_count]
        actual_mean[:, arm] = abs(data['nuisance_drift'][:, indices, arm].mean(axis=1))
    for name, budget in configurations:
        results[name] = infer_repair(data, target_calibration='empirical', target_alpha=.02,
            bounded=True, variance_failure=.005, bias_failure=.005,
            source_mode='dispersion', dispersion_calibration='bounded_exponential',
            source_bias_budget=budget['source_mean_budget'], target_bias_budget=budget['target_budget'])
        diagnostics[name] = dict(source_budget=budget['source_mean_budget'].mean(axis=0).tolist(),
            target_budget=budget['target_budget'].mean(axis=0).tolist(),
            mean_drift=actual_mean.mean(axis=0).tolist(),
            certificate_failure=float(np.mean(np.any(actual_mean > budget['source_mean_budget'] + 1e-12, axis=1)
                | np.any(data['oracle_target_budget'] > budget['target_budget'] + 1e-12, axis=1))))
    np.savez_compressed(str(Path(output)/'diagnostics'/('cell_%03d.npz'%setting['cell_id'])),
        methods=np.array(list(results)), source_budget=np.stack([b['source_mean_budget'] for _,b in configurations],axis=1),
        target_budget=np.stack([b['target_budget'] for _,b in configurations],axis=1),
        actual_mean_drift=actual_mean, target_reference_interval=reference)
    (Path(output)/'diagnostics'/('cell_%03d.json'%setting['cell_id'])).write_text(json.dumps(diagnostics,indent=2)+'\n')
    return save_case(output, setting, data, results, time.time()-started)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output')
    parser.add_argument('--repeats', type=int, default=300)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--reuse-source-design', action='store_true')
    parser.add_argument('--confirmation', action='store_true')
    parser.add_argument('--large-k', action='store_true')
    args = parser.parse_args()
    output = Path(args.output)
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/') or output.exists():
        raise ValueError('Use a new FACE-HD scratch output directory')
    if args.repeats < 1 or not 1 <= args.workers <= 4:
        raise ValueError('Positive repeats and one to four workers required')
    output.mkdir(parents=True)
    for name in ['cells','diagnostics','envelopes']:
        (output/name).mkdir()
    plan = make_plan()
    if args.confirmation:
        plan = [dict(cell_id=i+1, panel='fitted_confirmation', source_count=count,
            pattern=pattern, shared_scale=.3, n_per_site=1000,
            source_training_size=500, source_evaluation_size=250, target_evaluation_size=500,
            pilot_size=250, valid_fraction=.75, seed_offset=11130000,
            certified_reference='full_target_conditional')
            for i, (count, pattern) in enumerate([(128, 'weak_boundary'), (128, 'strong'),
                                                 (512, 'weak_boundary'), (512, 'strong')])]
    if args.large_k:
        if args.confirmation:
            raise ValueError('Choose confirmation or large-K exploration')
        plan = [dict(cell_id=i+1, panel='fitted_large_k_exploration', source_count=2048,
            pattern=pattern, shared_scale=.3, n_per_site=1000,
            source_training_size=500, source_evaluation_size=250, target_evaluation_size=500,
            pilot_size=250, valid_fraction=.75, seed_offset=12130000,
            certified_reference='full_target_conditional')
            for i, pattern in enumerate(['all_valid', 'weak_boundary', 'strong'])]
    if args.reuse_source_design:
        for setting in plan:
            setting.update(reuse_source_design=True, source_evaluation_size=500,
                           source_design_sample_size=500, source_design_evaluation_overlap=500,
                           pilot_size=0)
    if args.smoke:
        plan = plan[:1]
    configuration = dict(settings=plan,repeats=args.repeats,
        probability_lower_bounds=[.02,.04,.06],grid_size=2048,
        scope='Fixed pre-sampling valid sets; four-stratum design-only weights; MGF variants assume source joint-cell lower bounds. No production high-dimensional claim.',
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')})
    (output/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    envelopes={}
    for key,pmin in [('002',.02),('004',.04),('006',.06)]:
        envelopes[key]=binomial_ratio_mgf(500,4,pmin)
        np.savez_compressed(str(output/'envelopes'/('p'+key+'.npz')),**envelopes[key])
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for rows in pool.imap_unordered(run_cell,[(s,args.repeats,str(output),envelopes) for s in plan]):
            records.extend(rows)
            print('Completed cell',rows[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(sorted(records,key=lambda r:(r['cell_id'],r['method'])))
    (output/'COMPLETE').write_text('All paired samples and outcomes retained.\n')


if __name__=='__main__':
    main()

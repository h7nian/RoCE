#!/usr/bin/env python3
"""Export a complete paired p100 or p200 MC200 TATE panel without refitting."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

METHODS = {
    'one_round_crossfit_ate_joint_tate': 'RoCE (joint TATE)',
    'target_only_ate': 'Target-only',
    'target_anchor_ate': 'Calibrated target-only',
    'sample_size_ate': 'SS',
    'inverse_variance_ate': 'IVW',
    'federated_dr_ate': 'Federated-DR',
    'pooled_dr_ate': 'Pooled-DR'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(report, output):
    checks = json.loads((report/'pairing_checks.json').read_text())
    dimension = checks['dimension']
    if (checks['paired_repeats'] != 12000 or dimension not in (100, 200) or
            checks['crossfit_layers'] != 2 or checks['data_hash_matches'] != 12000 or
            checks['target_reference_matches'] != 12000):
        raise ValueError('Require a complete verified p100 or p200 Two-layer panel')
    with (report/'metrics.csv').open() as stream:
        metrics = [row for row in csv.DictReader(stream) if row['method'] in METHODS]
    with (report/'paired_estimates.csv').open() as stream:
        groups = defaultdict(dict)
        for row in csv.DictReader(stream):
            if row['method'] not in METHODS:
                continue
            if int(row['p']) != dimension:
                raise ValueError('Paired record dimension differs from the verified panel')
            key = row['config'], int(row['K']), float(row['rho']), row['method']
            seed = int(row['sim_id'])
            if seed in groups[key]:
                raise ValueError('Duplicate method/seed record')
            groups[key][seed] = row
    expected = {(config, K, rho, method) for config in ('C1','C2','C3')
                for K in (2,4,6,8) for rho in (0,.5,1,1.5,2) for method in METHODS}
    if set(groups) != expected or len(metrics) != len(expected):
        raise ValueError('Incomplete scenario/source/departure panel')
    exported = []
    for row in metrics:
        if int(row['p']) != dimension:
            raise ValueError('Metric dimension differs from the verified panel')
        key = row['config'], int(row['K']), float(row['rho']), row['method']
        records = groups[key]
        if set(records) != set(range(1,201)) or int(row['repeats']) != 200:
            raise ValueError('Every cell must retain all200 prescribed seeds')
        errors = [float(records[seed]['estimate'])-float(records[seed]['truth']) for seed in range(1,201)]
        squared = [value*value for value in errors]
        rmse = math.sqrt(statistics.mean(squared))
        if abs(rmse-float(row['rmse'])) > 1e-12:
            raise ValueError('RMSE does not reproduce the paired report')
        for seed in range(1,201):
            hashes = {groups[key[:3]+(method,)][seed]['data_sha256'] for method in METHODS}
            if len(hashes) != 1:
                raise ValueError('Figure methods use different generated data')
        mcse = statistics.stdev(squared)/math.sqrt(200)/(2*rmse)
        exported.append(dict(row, method_label=METHODS[row['method']], rmse_mcse=mcse,
                             rmse_lower=max(0,rmse-1.96*mcse),rmse_upper=rmse+1.96*mcse))
    output.mkdir(parents=True,exist_ok=False)
    with (output/'simulation_joint_tate_mc200.csv').open('x',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(exported[0]))
        writer.writeheader();writer.writerows(exported)
    # Keep paper-side provenance anonymous: detailed absolute run paths stay in scratch.
    provenance = dict(dgp_version='bounded_joint_v3',p=dimension,source_counts=[2,4,6,8],
        n_per_site=1000,scenarios=['C1','C2','C3'],rho=[0,.5,1,1.5,2],repeats_per_cell=200,
        crossfit_layers=2,communication='one_round',aggregation='joint_tate',cutoff=2,
        nuisance_lambdas=100,solver='proximal_newton',nuisance_tolerance=1e-10,
        predictor_radius=12,covariate_shift=1,source_treatment_scale=1,
        methods=METHODS,paired_data_settings=12000,metric_rows=len(exported),
        report_hashes={name:digest(report/name) for name in ['pairing_checks.json','metrics.csv','paired_estimates.csv']},
        monte_carlo_uncertainty='Wilson coverage intervals; delta-method MCSE for RMSE. Same generated data across methods. No failed seed removed.')
    (output/'simulation_joint_tate_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Exported',len(exported),'metric rows from12000 paired data settings.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    export(args.report,args.output)

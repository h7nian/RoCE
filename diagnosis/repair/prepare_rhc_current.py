#!/usr/bin/env python3
"""Freeze the current RHC application and predefined sensitivities for Slurm."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import EXCLUDED_NODES

SOURCE_RULE_VARIANTS = {
    'min': None,
    'final_weight_1se': dict(weight='1se'),
    'final_outcome_1se': dict(outcome='1se'),
    'initial_weight_1se': dict(initial_weight='1se'),
    'source_all_1se': dict(initial_weight='1se', weight='1se', outcome='1se'),
    'global_1se': None,
}
RHC_TARGET_SITES = ('Private', 'Medicare', 'Private & Medicare', 'Medicaid')


def source_rule_grid(fold_seeds=None, variants=None):
    """Resolve an explicit paired grid without changing the original 24-fit design."""
    seeds = [0, 101, 202, 303] if fold_seeds is None else list(fold_seeds)
    names = list(SOURCE_RULE_VARIANTS) if variants is None else list(variants)
    if (not seeds or any(type(seed) is not int or not 0 <= seed <= 2147483647 for seed in seeds)
            or len(set(seeds)) != len(seeds)):
        raise ValueError('Use distinct nonnegative integer fold seeds within the R integer range')
    if not names or len(set(names)) != len(names) or not set(names).issubset(SOURCE_RULE_VARIANTS):
        raise ValueError('Use distinct named source-rule variants')
    return [(seed, name, SOURCE_RULE_VARIANTS[name]) for seed in seeds for name in names]


def prepare(output, check_root, interface_checks, profiles=None, covariate_profile='historical',
            memory='24G', time_limit='06:00:00', study_profile='standard',
            fold_seeds=None, source_variants=None, target_site='Private', source_radii=None):
    if target_site not in RHC_TARGET_SITES:
        raise ValueError('Select one of the four retained RHC insurance sites as target')
    if study_profile not in ('standard', 'source_regularization', 'source_truncation'):
        raise ValueError('Unknown RHC study profile')
    if study_profile == 'standard' and (fold_seeds is not None or source_variants is not None):
        raise ValueError('Custom folds require a source-rule or source-truncation study')
    if source_radii is not None and study_profile != 'source_truncation':
        raise ValueError('Source radii require the source_truncation profile')
    if study_profile == 'source_truncation' and source_variants is not None:
        raise ValueError('Source truncation keeps the min rule at every stage')
    rule_grid = source_rule_grid(fold_seeds, source_variants) if study_profile == 'source_regularization' else None
    if study_profile == 'source_truncation':
        source_radii = [2, 3, 4, 5] if source_radii is None else list(source_radii)
        if not source_radii or len(set(source_radii)) != len(source_radii) or any(type(r) is not int or r not in (2, 3, 4, 5) for r in source_radii):
            raise ValueError('Select distinct prespecified source radii from 2,3,4,5')
        truncation_seeds = list(dict.fromkeys(seed for seed, _, _ in source_rule_grid(fold_seeds, ['min'])))
    output, check_root, interface_checks = map(pilot.scratch_path, (output, check_root, interface_checks))
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('RHC studies must remain in the real_data/rhc scratch hierarchy')
    markers = ('TESTS_PASSED', 'CROSSFIT_LAYER_CHECKS_PASSED', 'CV_CERTIFICATE_CHECKS_PASSED')
    for name in markers:
        if not (check_root / name).is_file():
            raise ValueError('Missing frozen-library check: ' + name)
    if study_profile == 'source_regularization' and not (check_root/'SOURCE_LAMBDA_RULES_CHECKS_PASSED').is_file():
        raise ValueError('Source-stage rule comparisons require their regression gate')
    if not (interface_checks / 'INTERFACE_CHECKS_PASSED').is_file():
        raise ValueError('The current RHC interface must pass its regression checks')
    if study_profile == 'source_truncation' and not (interface_checks/'SOURCE_TRUNCATION_CHECKS_PASSED').is_file():
        raise ValueError('Source truncation requires its target-isolation check')
    for name, digest in json.loads((interface_checks / 'checked_files.json').read_text()).items():
        if pilot.digest(Path(name)) != digest:
            raise ValueError('The RHC interface changed after its tests: ' + name)
    source = Path(__file__).resolve().parent
    repository = source.parent.parent
    library = check_root / 'Rlib'
    data_path = library / 'RoCE/extdata/rhc.csv'
    data_hash = pilot.digest(data_path)
    if data_hash != '9ef4ab578be4b40ad5d97d3a7e08ffdc1f9f76aeeefee51b4996e4221556f8e8':
        raise ValueError('RHC data differ from the reviewed public cohort')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('workflow', 'logs', 'tasks', 'submissions'):
        (output / name).mkdir()
    for name in ('submit_repeat_pilot.py', 'run_repeat_pilot.sh', 'run_rhc_current.R',
                 'layer_comparison_summary.R', 'rhc_stage_signatures.R', 'prepare_rhc_current.py'):
        shutil.copy2(str(source / name), str(output / 'workflow' / name))
    shutil.copy2(str(repository / 'R/real_data_rhc.R'), str(output / 'workflow/real_data_rhc.R'))
    tasks = []

    def add(name, role='method', baseline='', radius=5, nuisance_rule='min',
            source_program='calibrated', target_program='hou_calibrated', fold_seed=0,
            source_rules=None):
        tasks.append(dict(task_id=len(tasks)+1, config='RHC_'+name, K=3, sim_id=42,
            role=role, baseline=baseline, radius=radius, nuisance_rule=nuisance_rule,
            source_program=source_program, target_program=target_program, fold_seed=fold_seed))
        if study_profile == 'source_regularization':
            rules = source_rules or {}
            tasks[-1].update(source_initial_weight_rule=rules.get('initial_weight',''),
                             source_weight_rule=rules.get('weight',''), source_outcome_rule=rules.get('outcome',''))
        if study_profile == 'source_truncation':
            tasks[-1]['target_radius'] = 5

    if study_profile == 'standard':
        add('primary')
        for method in pilot.BASELINE_METHODS:
            add(method, role='baseline', baseline=method)
        add('nuisance_1se', nuisance_rule='1se')
        add('radius12', radius=12)
        add('standard_source', source_program='standard')
        add('ordinary_target', target_program='lasso')
        for seed in (101, 202, 303):
            add('fold_'+str(seed), fold_seed=seed)
    elif study_profile == 'source_regularization':
        for seed, name, rules in rule_grid:
            add(name+'_fold'+str(seed), fold_seed=seed,
                nuisance_rule='1se' if name=='global_1se' else 'min', source_rules=rules)
    else:
        for seed in truncation_seeds:
            for radius in source_radii:
                add('source_radius{}_fold{}'.format(radius, seed), radius=radius, fold_seed=seed)
    if profiles is not None:
        if not profiles or len(set(profiles)) != len(profiles) or not set(profiles).issubset(row['config'] for row in tasks):
            raise ValueError('Select distinct existing RHC profile names')
        tasks = [row for row in tasks if row['config'] in profiles]
    with (output / 'manifest.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(tasks[0]))
        writer.writeheader()
        writer.writerows(tasks)
    configuration = dict(scope='RHC current-method application and prespecified numerical/split sensitivities; no known causal truth',
        library=str(library), check_root=str(check_root), interface_checks=str(interface_checks),
        data_sha256=data_hash, worker_script='run_rhc_current.R', cpus=4, memory=memory,
        time_limit=time_limit, checkpoint=True, max_restarts=20, account='hou00123',
        partitions=pilot.DEFAULT_PARTITIONS, exclude_nodes=EXCLUDED_NODES,
        crossfit_layers=2, baseline_methods=list(pilot.BASELINE_METHODS),
        selected_profiles=[row['config'] for row in tasks],
        covariate_profile=covariate_profile,
        study_profile=study_profile,
        n_sites=4, source_count=3, target_site=target_site, n_folds=10, nlambda=100,
        n_bootstrap=5000, aggregation_cutoff=2, preprocessing='Historical pooled covariate-only imputation and standardization; shared OR/weight feature map; not fold-local preprocessing',
        communication_scope='Retrospective federated estimation emulation after pooled covariate preprocessing; not an end-to-end privacy implementation',
        inferential_scope='Observational treatment contrast conditional on identification/transport assumptions; split seeds do not define independent cohorts or coverage')
    pilot.write_json(output / 'configuration.json', configuration)
    design_name = {'standard': 'real_data_current_design.md',
                  'source_regularization': 'rhc_source_regularization_design.md',
                  'source_truncation': 'rhc_source_truncation_design.md'}[study_profile]
    shutil.copy2(str(source / design_name), str(output / 'README.md'))
    if profiles is not None:
        with (output / 'README.md').open('a') as stream:
            stream.write('\n## This submission subset\n\nSelected profiles: '+', '.join(profiles)+'. Original task IDs are retained.\n')
    with (output / 'README.md').open('a') as stream:
        stream.write('\nCovariate profile: `'+covariate_profile+'`. Memory: '+memory+'; wall time: '+time_limit+'.\n')
        stream.write('\nTarget site: `'+target_site+'`; the other three retained insurance groups are sources. '
                     'This configured target takes precedence over the historical Private-target design above.\n')
        if study_profile == 'source_regularization':
            stream.write('\nThis frozen manifest contains {} fits; fold seeds: {}; variants: {}. '
                         'It takes precedence over the original discovery grid described above.\n'.format(
                             len(tasks), ', '.join(map(str, dict.fromkeys(row['fold_seed'] for row in tasks))),
                             ', '.join(dict.fromkeys(name for _, name, _ in rule_grid))))
    if covariate_profile != 'historical':
        shutil.copy2(str(source / 'rhc_covariate_diagnosis.md'), str(output / 'covariate_design.md'))
    paths = [output / 'configuration.json', output / 'manifest.csv', output / 'README.md',
             interface_checks / 'INTERFACE_CHECKS_PASSED', interface_checks / 'checked_files.json']
    paths += [check_root / name for name in markers]
    if covariate_profile != 'historical':
        paths.append(output / 'covariate_design.md')
    if study_profile == 'source_regularization':
        paths.append(check_root/'SOURCE_LAMBDA_RULES_CHECKS_PASSED')
    if study_profile == 'source_truncation':
        paths.append(interface_checks/'SOURCE_TRUNCATION_CHECKS_PASSED')
    paths += sorted((output / 'workflow').iterdir())
    paths += sorted(path for path in (library / 'RoCE').rglob('*') if path.is_file())
    pilot.write_json(output / 'provenance.json', {str(path): pilot.digest(path) for path in paths})
    print('Prepared {} RHC jobs ({} method profiles and {} baselines).'.format(
        len(tasks), sum(row['role']=='method' for row in tasks), sum(row['role']=='baseline' for row in tasks)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--check-root', required=True, type=Path)
    parser.add_argument('--interface-checks', required=True, type=Path)
    parser.add_argument('--profiles', nargs='+', help='Submit only these named profiles, retaining original task IDs')
    parser.add_argument('--covariate-profile', choices=('historical', 'log_missing', 'grouped_log_missing'), default='historical')
    parser.add_argument('--memory', default='24G')
    parser.add_argument('--time-limit', default='06:00:00')
    parser.add_argument('--study-profile', choices=('standard','source_regularization','source_truncation'), default='standard')
    parser.add_argument('--fold-seeds', nargs='+', type=int,
                        help='Source-rule/truncation studies; 0 retains the historical default partition')
    parser.add_argument('--source-variants', nargs='+', choices=tuple(SOURCE_RULE_VARIANTS),
                        help='Source-rule study only; omitted selects all six original variants')
    parser.add_argument('--target-site', choices=RHC_TARGET_SITES, default='Private')
    parser.add_argument('--source-radii', nargs='+', type=int, choices=(2,3,4,5),
                        help='Source-truncation profile only; target radius remains 5')
    args = parser.parse_args()
    prepare(args.output, args.check_root, args.interface_checks, args.profiles,
            args.covariate_profile, args.memory, args.time_limit, args.study_profile,
            args.fold_seeds, args.source_variants, args.target_site, args.source_radii)

#!/usr/bin/env python3
"""Freeze targeted RHC case diagnostics with unchanged outer-fold memberships."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import EXCLUDED_NODES


def reference_directory(base, seed, radius):
    if seed == 0:
        study = base/'source_truncation_v1'
        profile = 'RHC_source_radius{}_fold0'.format(radius)
    elif radius == 5:
        study = base/'source_regularization_confirmation_v1'
        profile = 'RHC_min_fold'+str(seed)
    else:
        study = base/'source_truncation_confirmation_v1'
        profile = 'RHC_source_radius{}_fold{}'.format(radius, seed)
    with (study/'manifest.csv').open() as stream:
        rows = [r for r in csv.DictReader(stream) if r['config'] == profile]
    if len(rows) != 1:
        raise ValueError('Missing unique saved reference: '+profile)
    directory = study/'tasks'/rows[0]['task_id']
    if not (directory/'COMPLETE').is_file():
        raise ValueError('Reference is incomplete: '+profile)
    return directory


def prepare(output, checks):
    output, checks = pilot.scratch_path(output), pilot.scratch_path(checks)
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Use the RHC scratch hierarchy')
    if not (checks/'CHECKS_PASSED').is_file():
        raise ValueError('Case/fold construction checks must pass')
    for name, expected in json.loads((checks/'checked_files.json').read_text()).items():
        if pilot.digest(Path(name)) != expected:
            raise ValueError('Case-refit input code changed after checks: '+name)
    source = Path(__file__).resolve().parent
    base = output.parent
    library = Path('/scratch.global/zhan9381/FACE-HD/implementation/r12/source_lambda_rules_v2/Rlib')
    output.mkdir(parents=True, exist_ok=False)
    for directory in ('workflow', 'tasks', 'logs', 'submissions'):
        (output/directory).mkdir()
    for name in ('submit_repeat_pilot.py', 'run_repeat_pilot.sh', 'rhc_validation_helpers.R',
                 'run_rhc_case_refit.R', 'rhc_stage_signatures.R', 'layer_comparison_summary.R'):
        shutil.copy2(str(source/name), str(output/'workflow'/name))
    tasks = []
    # Published partition and worst observed partitions for radius2/3/5.
    cases = [(0, 's3', 423), (412, 's3', 423), (419, 't', 303)]
    for seed, site, record in [(s, 'none', 0) for s in (0, 412, 419)] + cases:
        for radius in (2, 3, 5):
            tasks.append(dict(task_id=len(tasks)+1, K=3, sim_id=seed,
                config='RHC_case_{}_{}_radius{}_fold{}'.format(site, record, radius, seed),
                fold_seed=seed, case_site=site, case_row=record, radius=radius,
                reference=str(reference_directory(base, seed, radius))))
    with (output/'manifest.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(tasks[0])); writer.writeheader(); writer.writerows(tasks)
    configuration = dict(library=str(library), worker_script='run_rhc_case_refit.R',
        cpus=4, memory='8G', time_limit='01:00:00', checkpoint=True, max_restarts=20,
        account='hou00123', partitions=pilot.DEFAULT_PARTITIONS, exclude_nodes=EXCLUDED_NODES,
        crossfit_layers=2, target='Private', target_radius=5, source_radii=[2,3,5],
        nuisance_rule='min', aggregation_cutoff=2,
        scope='Full nuisance/CV/eta refit after one diagnostic removal; original preprocessing and remaining outer folds fixed')
    pilot.write_json(output/'configuration.json', configuration)
    (output/'README.md').write_text('''# RHC complete-refit case sensitivity

The retained cohort,Private target,min rule,Two-layer/one-round,joint TATE and
target radius5 stay fixed. Compare source radii2,3,5. Preserve every remaining
patient's outer-fold assignment andthe original common preprocessing; retune
andrefit all nuisance models andeta. Nuisance CV is rerun by the original rule,
so this measures the full fitting program's response,not a fixed-model influence.

Three prespecified case/partition pairs: published partition0/Medicaid423;
partition412/Medicaid423 (largest recorded influence under radius3 andradius5);
partition419/target303 (largest recorded influence under radius2). Run each
radius onthe full cohort andafter each removal:9 controls plus9 perturbations.
First submit only the9 full-cohort controls andrequire exact saved-fit parity.

Report every case,including adverse changes,as a diagnostic. Do not delete a
patient from the primary analysis or select a favorable leave-out result.
These checks do not identify causal bias or supply a new clinical estimate.
All original data/results remain unchanged.
''')
    paths = [output/'manifest.csv', output/'configuration.json', output/'README.md', checks/'CHECKS_PASSED', checks/'checked_files.json']
    paths += list((output/'workflow').iterdir())
    paths += [p for p in (library/'RoCE').rglob('*') if p.is_file()]
    for task in tasks:
        paths += [Path(task['reference'])/'COMPLETE', Path(task['reference'])/'checkpoints/analysis.rds']
    pilot.write_json(output/'provenance.json', {str(p):pilot.digest(p) for p in paths})
    print('Prepared18 controlled case fits; controls are tasks1–9.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--checks', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.output, args.checks)

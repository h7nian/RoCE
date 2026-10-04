#!/usr/bin/env python3
"""Freeze failed-repeat recovery and four successful-path controls for RHC."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot


def prepare(output, parent, library):
    output, parent, library = map(pilot.scratch_path, (output, parent, library))
    configuration = pilot.verify(parent)
    checks = library.parent
    for name in ('DENSITY_CV_RETRY_CHECKS_PASSED', 'RHC_CV_FIXTURES_PASSED'):
        if not (checks/name).is_file():
            raise ValueError('Recovery gate missing: '+name)
    with (parent/'manifest.csv').open() as stream:
        manifest = list(csv.DictReader(stream))
    jobs = [json.loads((parent/'submissions'/(row['task_id']+'.json')).read_text())['stdout'].strip().split(';')[0]
            for row in manifest]
    live = subprocess.check_output(['squeue', '-j', ','.join(jobs), '-h', '-o', '%i'], universal_newlines=True)
    if live.strip():
        raise ValueError('Parent scientific jobs are still live; do not freeze an incomplete recovery list')
    failed = []
    for row in manifest:
        directory = parent/'tasks'/row['task_id']
        if (directory/'COMPLETE').is_file():
            continue
        failure = directory/'analysis_FAILED.txt'
        if not failure.is_file() or 'no lambda converged with a finite validation loss across every CV fold' not in failure.read_text():
            raise ValueError('Unclassified failure requires inspection: '+row['task_id'])
        failed.append(row['task_id'])
    if not failed:
        raise ValueError('No failed repeats require recovery')
    controls = [row['task_id'] for row in manifest[:4]]
    if set(failed).intersection(controls):
        raise ValueError('The first four successful control repeats are required')
    selected = [row for row in manifest if row['task_id'] in failed + controls]
    source = Path(__file__).resolve().parent
    output.mkdir(parents=True, exist_ok=False)
    for name in ('workflow', 'tasks', 'submissions', 'logs'):
        (output/name).mkdir()
    for name in ('submit_repeat_pilot.py', 'run_repeat_pilot.sh', 'rhc_validation_helpers.R',
                 'rhc_model_validation.R', 'rhc_stage_signatures.R', 'run_rhc_model_validation.R',
                 'review_rhc_model_validation.R', 'review_rhc_model_recovery.R'):
        shutil.copy2(str(source/name), str(output/'workflow'/name))
    with (output/'manifest.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest[0])); writer.writeheader(); writer.writerows(selected)
    configuration.update(library=str(library), parent_study=str(parent),
        recovery_task_ids=failed, control_task_ids=controls,
        scope='Fresh refits of all failed original repeats plus successful-path controls; same populations, seeds and fitting controls')
    pilot.write_json(output/'configuration.json', configuration)
    (output/'README.md').write_text('''# RHC model validation: failed-repeat recovery

Keep the original four population laws, 20 repetitions per law, exact seeds,
site sizes, feature basis, min CV rule, two-layer/one-round joint TATE and radii
2/3/5. The shared density-CV retry changes only entirely failed CV paths.

Refit all failed repeats from scratch. Also refit the four originally successful
first repeats as controls. Do not migrate caches across installed-code versions.
Require successful controls to reproduce the saved coefficients, aggregation
weights, point estimates and SEs within 1e-10, and require fresh real-RHC radius
3/5 parity. Only then combine repaired repeats with the unchanged successful
parent results. The combined panel contains each original scenario/seed exactly
once. Report all original failures and their disposition; retain both versions.

No observed RHC outcome, target group, patient inclusion, population template or
replication seed is changed. Twenty repeats remain a model-based pilot, not proof
of real RHC causal bias or precise 95% coverage.
''')
    paths = [output/'manifest.csv', output/'configuration.json', output/'README.md',
             parent/'manifest.csv', parent/'configuration.json', parent/'provenance.json',
             checks/'DENSITY_CV_RETRY_CHECKS_PASSED', checks/'RHC_CV_FIXTURES_PASSED', Path(configuration['template'])]
    paths += list((output/'workflow').iterdir())
    paths += [p for p in (library/'RoCE').rglob('*') if p.is_file()]
    pilot.write_json(output/'provenance.json', {str(p):pilot.digest(p) for p in paths})
    print('Prepared {} recovered repeats and {} successful controls.'.format(len(failed), len(controls)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--parent', required=True, type=Path)
    parser.add_argument('--library', required=True, type=Path)
    args = parser.parse_args()
    prepare(args.output, args.parent, args.library)

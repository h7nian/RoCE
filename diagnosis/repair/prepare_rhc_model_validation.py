#!/usr/bin/env python3
"""Freeze paired known-truth RHC model diagnostics, one replicate per Slurm job."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import EXCLUDED_NODES


def prepare(output, template, checks, repeats, library, seed_start=1, scenarios=None):
    output, template, checks = map(pilot.scratch_path, (output, template, checks))
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Use the RHC scratch hierarchy')
    if (repeats < 2 or seed_start < 1 or seed_start + repeats - 1 > 2146544647 or
            not all((path/'CHECKS_PASSED').is_file() for path in (template, checks))):
        raise ValueError('Require at least two repeats and a checked population template')
    source = Path(__file__).resolve().parent
    scenarios = scenarios or ['O_case_mix', 'W_overlap', 'W_tail', 'W_departure']
    with (template/'population_truth.csv').open() as stream:
        available = {row['scenario'] for row in csv.DictReader(stream)}
    if len(set(scenarios)) != len(scenarios) or not set(scenarios).issubset(available):
        raise ValueError('Scenarios must be distinct and present in the checked template')
    library = pilot.scratch_path(library)
    if not (library/'RoCE'/'DESCRIPTION').is_file():
        raise ValueError('Installed library is missing')
    for directory in (template, checks):
        checked = json.loads((directory/'checked_files.json').read_text())
        for name, expected in checked.items():
            if pilot.digest(Path(name)) != expected:
                raise ValueError('Population code/input changed after checking: '+name)
    output.mkdir(parents=True, exist_ok=False)
    for name in ('workflow', 'tasks', 'submissions', 'logs'):
        (output/name).mkdir()
    for name in ('submit_repeat_pilot.py', 'run_repeat_pilot.sh', 'rhc_validation_helpers.R',
                 'rhc_model_validation.R', 'rhc_stage_signatures.R', 'run_rhc_model_validation.R',
                 'review_rhc_model_validation.R'):
        shutil.copy2(str(source/name), str(output/'workflow'/name))
    shutil.copy2(str(source/'rhc_model_validation_design.md'), str(output/'workflow/background_design.md'))
    design_link = 'workflow/background_design.md'
    if 'O_moderate_mix' in scenarios:
        shutil.copy2(str(source/'rhc_or_branch_design.md'), str(output/'workflow/or_branch_design.md'))
        design_link = 'workflow/or_branch_design.md'
    (output/'README.md').write_text('''# RHC model validation: frozen campaign

Selected population laws: {scenarios}.
Repeat IDs: {start}–{end} inclusive; {repeats} repeats per selected law.
Each repeat pairs source radii2,3,5 andretains target radius5.
Target-only denotes the calibrated anchor. Every original result is retained.

The selected [population design]({design}) andthe configuration/manifest define
this campaign. Background-family descriptions do not override the selected
scenarios orrepeat range. Results are model-based diagnostics,not observed RHC
bias estimates. A continuation range retains earlier prescribed repeats.
'''.format(scenarios=', '.join(scenarios),start=seed_start,end=seed_start+repeats-1,
           repeats=repeats,design=design_link))
    rows = []
    for repeat in range(seed_start, seed_start+repeats):
        for scenario in scenarios:
            rows.append(dict(task_id=len(rows)+1, config=scenario, K=3, sim_id=repeat,
                population_seed=929000+repeat, fit_seed=939000+repeat))
    with (output/'manifest.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    configuration = dict(library=str(library), template=str(template/'template.rds'),
        source_radii=[2, 3, 5], target_radius=5, repeats=repeats, seed_start=seed_start, scenarios=scenarios,
        worker_script='run_rhc_model_validation.R', cpus=4, memory='8G', time_limit='01:00:00',
        checkpoint=True, max_restarts=20, account='hou00123', partitions=pilot.DEFAULT_PARTITIONS,
        exclude_nodes=EXCLUDED_NODES, crossfit_layers=2,
        scope='Frozen model-based RHC truncation sensitivity; not actual RHC bias or formal MC coverage evidence')
    pilot.write_json(output/'configuration.json', configuration)
    paths = [output/'manifest.csv', output/'configuration.json', output/'README.md']
    paths += list((output/'workflow').iterdir()) + list(template.iterdir()) + list(checks.iterdir())
    paths += [p for p in (library/'RoCE').rglob('*') if p.is_file()]
    pilot.write_json(output/'provenance.json', {str(p):pilot.digest(p) for p in paths if p.is_file()})
    print('Prepared {} paired replicates for {} checked population laws.'.format(len(rows), len(scenarios)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--template', required=True, type=Path)
    parser.add_argument('--checks', required=True, type=Path)
    parser.add_argument('--repeats', default=20, type=int)
    parser.add_argument('--seed-start', default=1, type=int)
    parser.add_argument('--scenarios', nargs='+')
    parser.add_argument('--library', type=Path,
        default=Path('/scratch.global/zhan9381/FACE-HD/implementation/r12/source_lambda_rules_v2/Rlib'))
    args = parser.parse_args()
    prepare(args.output, args.template, args.checks, args.repeats, args.library, args.seed_start, args.scenarios)

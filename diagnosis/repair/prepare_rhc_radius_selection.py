#!/usr/bin/env python3
"""Select a saved Private/min radius and retain the matched baseline results."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def prepare(base, output, radius):
    base, output = base.resolve(), output.resolve()
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Use the RHC scratch hierarchy')
    previous = base/'report_selection_private_min_v1'
    candidates = [r for r in read_rows(base/'source_truncation_v1/manifest.csv')
                  if r['config'] == 'RHC_source_radius{}_fold0'.format(radius)]
    if len(candidates) != 1:
        raise ValueError('A unique original-partition fit is required')
    selected = base/'source_truncation_v1/tasks'/candidates[0]['task_id']
    if not (selected/'COMPLETE').is_file():
        raise ValueError('Selected fit is incomplete')
    reference = read_rows(previous/'methods.csv')
    replacements = {r['method']:r for r in read_rows(selected/'methods.csv')}
    for row in reference:
        if row['method'] in ('Target-only', 'Calibrated target-only'):
            if any(abs(float(row[k])-float(replacements[row['method']][k])) > 1e-12
                   for k in ('estimate', 'se', 'ci_lower', 'ci_upper')):
                raise ValueError('The target reference changed')
    output.mkdir(parents=True, exist_ok=False)
    with (output/'methods.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(reference[0])); writer.writeheader()
        writer.writerows([replacements.get(r['method'], r) for r in reference])
    for name in ('arm_summary.csv', 'sources.csv'):
        shutil.copy2(str(selected/name), str(output/name))
    configuration = json.loads((previous/'configuration.json').read_text())
    configuration.update(primary_fit=str(selected), source_radius=radius, target_radius=5,
        asset_prefix='rhc_private_min_radius'+str(radius),
        reporting_scope='Private/min, original partition, matched baselines, explicitly selected source radius',
        internal_only=['source_all_1se', 'global_1se', 'rotated_target_populations', 'unselected_covariate_panels'])
    (output/'configuration.json').write_text(json.dumps(configuration, indent=2)+'\n')
    paths = list(output.iterdir()) + [selected/'COMPLETE', selected/'checkpoints/analysis.rds',
        selected/'methods.csv', selected/'sources.csv', selected/'arm_summary.csv']
    provenance = json.loads((previous/'provenance.json').read_text())
    provenance.update({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--base', type=Path, default=Path('/scratch.global/zhan9381/FACE-HD/real_data/rhc'))
    parser.add_argument('--source-radius', type=int, choices=(2, 3, 4, 5), required=True)
    args = parser.parse_args()
    prepare(args.base, args.output, args.source_radius)

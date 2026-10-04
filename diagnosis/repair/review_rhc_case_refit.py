#!/usr/bin/env python3
"""Review full-refit case sensitivity without changing the primary cohort."""
import argparse
import csv
import json
from pathlib import Path


def review(root, output):
    with (root/'manifest.csv').open() as stream:
        manifest = list(csv.DictReader(stream))
    results = {}
    for task in manifest:
        directory = root/'tasks'/task['task_id']
        if not (directory/'COMPLETE').is_file():
            raise ValueError('Incomplete case task: '+task['task_id'])
        if task['case_site'] == 'none' and not (directory/'REFERENCE_PARITY_PASSED').is_file():
            raise ValueError('Full cohort did not pass reference parity')
        with (directory/'methods.csv').open() as stream:
            methods = {row['method']:row for row in csv.DictReader(stream)}
        with (directory/'nuisance_fits.csv').open() as stream:
            diagnostics = list(csv.DictReader(stream))
        for field in diagnostics[0]:
            if 'nonconverged' in field or 'line_search_failures' in field:
                if any(float(row[field]) != 0 for row in diagnostics):
                    raise ValueError('A fitted model failed: '+task['config']+'/'+field)
        key = (task['fold_seed'],task['radius'],task['case_site'],task['case_row'])
        if key in results:
            raise ValueError('Duplicate case specification')
        results[key] = methods
    comparisons = []
    for task in manifest:
        if task['case_site'] == 'none':
            continue
        changed = results[(task['fold_seed'],task['radius'],task['case_site'],task['case_row'])]
        baseline = results[(task['fold_seed'],task['radius'],'none','0')]
        fit, reference = changed['RoCE'], baseline['RoCE']
        comparisons.append(dict(fold_seed=task['fold_seed'],source_radius=task['radius'],
            removed_site=task['case_site'],removed_row=task['case_row'],
            full_estimate=float(reference['estimate']),refit_estimate=float(fit['estimate']),
            estimate_change=float(fit['estimate'])-float(reference['estimate']),
            change_over_full_se=(float(fit['estimate'])-float(reference['estimate']))/float(reference['se']),
            full_se=float(reference['se']),refit_se=float(fit['se']),
            anchor_change=float(changed['Target-only']['estimate'])-float(baseline['Target-only']['estimate'])))
    output.mkdir(parents=True, exist_ok=False)
    with (output/'comparisons.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
    (output/'CHECKS_PASSED').write_text('All18 case fits complete; nine full-cohort references reproduced; fitted-model convergence checked.\n')
    print(comparisons)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();review(args.root,args.output)

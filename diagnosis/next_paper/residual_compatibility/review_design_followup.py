"""Recompute paired uncertainty and create the v21 follow-up report assets."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta


def binomial_interval(successes, repeats):
    return (0. if successes == 0 else float(beta.ppf(.025, successes, repeats-successes+1)),
            1. if successes == repeats else float(beta.ppf(.975, successes+1, repeats-successes)))


def summarize_cell(path, cohort):
    record = json.loads(path.read_text())
    setting, names = record['setting'], record['methods']
    data = np.load(str(path.with_suffix('.npz')))
    intervals = data['interval']
    lengths = intervals[:, :, 1]-intervals[:, :, 0]
    coverage = (intervals[:, :, 0] <= .1) & (.1 <= intervals[:, :, 1])
    conditional = data['conditional']
    conditional_coverage = (conditional[:, :, 0] <= data['conditional_truth'][:, None]) & (data['conditional_truth'][:, None] <= conditional[:, :, 1])
    n = len(lengths)
    rows = []
    for j, method in enumerate(names):
        lower, upper = binomial_interval(int(coverage[:, j].sum()), n)
        # Independently check the stored summary against the saved interval endpoints.
        stored = next(m for m in record['metrics'] if m['method'] == method)
        np.testing.assert_allclose([stored['coverage'], stored['length'], stored['conditional_coverage']],
            [coverage[:, j].mean(), lengths[:, j].mean(), conditional_coverage[:, j].mean()], rtol=0, atol=1e-12)
        rows.append(dict(setting, cohort=cohort, method=method, repeats=n,
            length=float(lengths[:, j].mean()), length_mcse=float(lengths[:, j].std(ddof=1)/np.sqrt(n)),
            coverage=float(coverage[:, j].mean()), coverage_lower=lower, coverage_upper=upper,
            conditional_coverage=float(conditional_coverage[:, j].mean())))
    pairs = []
    if 'source_bernstein' in names:
        comparisons = [('source_bernstein', reference) for reference in ['target_normal', 'target_bernstein']]
    else:
        comparisons = [(method, reference) for method in names if method.startswith('remove_')
                       for reference in ['all_required_'+method[len('remove_'):], 'target_protected'] if reference in names]
    for method, reference in comparisons:
        j, k = names.index(method), names.index(reference)
        difference = lengths[:, j]-lengths[:, k]
        coverage_difference = coverage[:, j].astype(float)-coverage[:, k]
        pairs.append(dict(setting, cohort=cohort, method=method, reference=reference,
            repeats=n, length_difference=float(difference.mean()),
            paired_mcse=float(difference.std(ddof=1)/np.sqrt(n)),
            length_ratio=float(lengths[:, j].mean()/lengths[:, k].mean()),
            coverage_difference=float(coverage_difference.mean()),
            coverage_difference_mcse=float(coverage_difference.std(ddof=1)/np.sqrt(n))))
    return rows, pairs, record


def write_csv(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--previous', type=Path, required=True)
    args = parser.parse_args()
    report = args.root/'report'
    report.mkdir(exist_ok=True)
    rows, pairs, removal = [], [], []
    for folder, cohort in [(args.previous/'target', 'v20'), (args.root/'confirmation', 'independent'),
                            (args.root/'removal', 'removal')]:
        if not (folder/'COMPLETE').is_file():
            raise ValueError('Panel incomplete: '+str(folder))
        for path in sorted(folder.glob('cell_*.json')):
            current, comparisons, record = summarize_cell(path, cohort)
            rows.extend(current); pairs.extend(comparisons)
            if cohort == 'removal':
                audits = record['diagnostics']
                removal.append(dict(record['setting'],
                    any_failure=float(np.mean([r['rejected_count'] > 0 for r in audits])),
                    mean_retained=float(np.mean([r['retained_count'] for r in audits])),
                    eligibility=float(np.mean([r['eligible'] for r in audits])),
                    mean_max_weight_norm=float(np.mean([r['max_weight_square_sum'] for r in audits if r['max_weight_square_sum'] is not None]))))
    write_csv(report/'method_metrics.csv', rows)
    write_csv(report/'paired_comparisons.csv', pairs)
    write_csv(report/'removal_diagnostics.csv', removal)
    target_rows = []
    for cell in range(1, 5):
        chosen = [r for r in rows if r['cohort'] == 'independent' and r['cell_id'] == cell]
        by_name = {r['method']: r for r in chosen}
        s, t, protected = [by_name[name] for name in ['source_bernstein', 'target_normal', 'target_bernstein']]
        pair = next(r for r in pairs if r['cohort'] == 'independent' and r['cell_id'] == cell and r['reference'] == 'target_normal')
        target_rows.append('%d & %d & %.2f & %s & %.4f & %.4f & %.4f & %.3f/%.3f & %.4f (%.4f) \\\\' %
            (s['dimension'], s['source_count'], s['heterogeneity'], s['pattern'].replace('_', ' '),
             s['length'], protected['length'], t['length'], s['coverage'], t['coverage'],
             pair['length_difference'], pair['paired_mcse']))
    (report/'confirmation_table.tex').write_text('\\begin{tabular}{rrrlrrrrr}\n\\toprule\n'
        '$p$ & $K$ & $h$ & Pattern & Source & Protected target & Normal target & Coverage S/N & Difference (MCSE) \\\\\n\\midrule\n'
        +'\n'.join(target_rows)+'\n\\bottomrule\n\\end{tabular}\n')
    removal_rows = []
    for cell in [2, 4, 6, 8, 9, 10]:
        chosen = {r['method']: r for r in rows if r['cohort'] == 'removal' and r['cell_id'] == cell}
        audit = next(r for r in removal if r['cell_id'] == cell)
        removal_rows.append('%d & %d & %d & %.3f & %.3f & %.3f & %.4f & %.4f & %.4f \\\\' %
            (audit['source_size'], audit['dimension'], audit['source_count'], audit['singular_fraction'],
             audit['any_failure'], audit['eligibility'], chosen['all_required_linear']['length'],
             chosen['remove_linear']['length'], chosen['target_protected']['length']))
    (report/'removal_table.tex').write_text('\\begin{tabular}{rrrrrrrrr}\n\\toprule\n'
        '$n_s$ & $p$ & $K$ & Singular fraction & Any failure & Eligible & Require all & Remove & Target \\\\\n\\midrule\n'
        +'\n'.join(removal_rows)+'\n\\bottomrule\n\\end{tabular}\n')
    boundary = []
    for path in sorted((args.root/'model_boundary').glob('k*_eta*.json')):
        result = json.loads(path.read_text())
        data = np.load(str(path.with_suffix('.npz')))
        for j, row in enumerate(result['metrics']):
            interval = data['conditional_interval'][:, j]
            np.testing.assert_allclose(row['conditional_length'], np.diff(interval, axis=1).mean(), atol=1e-12)
            boundary.append(row)
    if len(boundary) != 16 or not (args.root/'model_boundary'/'COMPLETE').is_file():
        raise ValueError('Boundary panel incomplete')
    write_csv(report/'boundary_metrics.csv', boundary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.alpha': .2, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, eta in zip(axes, [0., .4]):
        for name, label, color in [('linear_basis_bounded_invalid', 'Bounded invalid, basis (1,x)', '#D55E00'),
                                   ('saturated_basis', 'Saturated basis (1,x,x²)', '#0072B2')]:
            values = sorted([r for r in boundary if r['method'] == name and r['nonlinearity'] == eta], key=lambda r:r['source_count'])
            ax.loglog([r['source_count'] for r in values], [r['source_length'] for r in values], '-o', color=color, label=label)
        ax.set_title('Nonlinearity amplitude '+str(eta))
        ax.set_xlabel('Number of sources K')
        ax.set_ylabel('Mean conditional source-band length')
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(str(report/'model_boundary.pdf')); fig.savefig(str(report/'model_boundary.png'), dpi=160)
    plt.close(fig)
    summary = dict(outcome_setting_repetitions=4400, confirmation_rows=[r for r in rows if r['cohort'] == 'independent'
        and r['method'] in ['source_bernstein', 'target_bernstein', 'target_normal']],
        key_pairs=[r for r in pairs if r['cohort'] == 'independent'], removal=removal,
        boundary=boundary, exact_lower_bounds=json.loads((args.root/'model_boundary'/'contamination_overlap.json').read_text()))
    (report/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('Reviewed', len(rows), 'method rows and', len(pairs), 'paired comparisons')


if __name__ == '__main__':
    main()

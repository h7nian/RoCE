"""Recompute continuous-balance metrics and produce the paired research report."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def write_table(path, columns, rows):
    path.write_text('\n'.join(['\\begin{tabular}{'+columns+'}\\toprule']+rows+
                             ['\\bottomrule', '\\end{tabular}'])+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root, rows, audits, paired = args.root, [], [], []
    report = root/'report'
    report.mkdir(exist_ok=True)
    for profile in ['main', 'boundary', 'stress']:
        if not (root/profile/'COMPLETE').exists():
            raise RuntimeError('Finish the prescribed profile before reporting: '+profile)
        for path in sorted((root/profile).glob('cell_*.json')):
            record = json.loads(path.read_text())
            arrays = np.load(str(path.with_suffix('.npz')))
            names = record['methods']
            for index, metric in enumerate(record['metrics']):
                interval = arrays['intervals'][:, index]
                covered = (interval[:, 0] <= .1) & (.1 <= interval[:, 1])
                length = np.diff(interval, axis=1)[:, 0]
                error = arrays['points'][:, index]-.1
                calculated = dict(coverage=covered.mean(), length=length.mean(),
                                  rmse=np.sqrt(np.mean(error**2)), bias=error.mean())
                for name, value in calculated.items():
                    if abs(value-metric[name]) > 1e-12:
                        raise AssertionError('Saved metric differs: '+str(path)+' '+name)
                count, total = int(covered.sum()), len(covered)
                metric.update(profile=profile, length_mcse=float(length.std(ddof=1)/np.sqrt(total)),
                    coverage_lower=float(beta.ppf(.025, count, total-count+1)) if count else 0.,
                    coverage_upper=float(beta.ppf(.975, count+1, total-count)) if count < total else 1.)
                rows.append(metric)
            for suffix in ['universal', 'declared_half']:
                target = names.index('target_'+suffix)
                base_length = np.diff(arrays['intervals'][:, target], axis=1)[:, 0]
                for mode in ['arms', 'adaptive']:
                    method = mode+'_'+suffix
                    index = names.index(method)
                    difference = np.diff(arrays['intervals'][:, index], axis=1)[:, 0]-base_length
                    mse_change = (arrays['points'][:, index]-.1)**2-(arrays['points'][:, target]-.1)**2
                    paired.append(dict(profile=profile, cell_id=record['setting']['cell_id'], method=method,
                        length_change=float(difference.mean()),
                        length_change_mcse=float(difference.std(ddof=1)/np.sqrt(len(difference))),
                        mse_change=float(mse_change.mean()),
                        mse_change_mcse=float(mse_change.std(ddof=1)/np.sqrt(len(mse_change)))))
            diagnostics = record['diagnostics']
            supported = [r for r in diagnostics if not r['source_failed']]
            audits.append(dict(profile=profile, **record['setting'], repeats=len(diagnostics),
                source_fallback=sum(r['source_failed'] for r in diagnostics),
                target_fallback=sum(r['target_failed'] for r in diagnostics),
                max_balance_error=max([r.get('maximum_balance_error', 0.) for r in diagnostics]),
                max_valid_drift=max([r.get('maximum_valid_drift', 0.) for r in diagnostics]),
                max_leverage=max([r.get('max_leverage', 1.) for r in diagnostics]),
                negative_weight_fraction=float(np.mean([r['negative_weight_fraction'] for r in supported])) if supported else None,
                dispersion_failure=sum(r.get('dispersion_failures', 0) for r in diagnostics),
                naive_conditional_coverage=float(np.mean([r['naive_conditional_coverage'] for r in diagnostics])),
                variance_ratio=float(np.mean([r['source_variance_estimate'] for r in supported]) /
                                     np.mean([r['source_variance_oracle'] for r in supported])) if supported else None,
                runtime_seconds=record['runtime_seconds']))
    for filename, values in [('metrics.csv', rows), ('paired_comparisons.csv', paired), ('design_audit.csv', audits)]:
        with (report/filename).open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0])); writer.writeheader(); writer.writerows(values)
    lookup = {(r['profile'], r['cell_id'], r['method']): r for r in rows}
    protected = [r for r in rows if r['method'] != 'target_normal' and r['pattern'] != 'nonlinear']
    exact = [r for r in audits if r['pattern'] != 'nonlinear']
    summary = dict(distinct_setting_repetitions=sum(r['repeats'] for r in audits),
        reported_settings=len(audits), method_setting_rows=len(rows),
        protected_coverage_range=[min(r['coverage'] for r in protected), max(r['coverage'] for r in protected)],
        maximum_linear_valid_drift=max(r['max_valid_drift'] for r in exact),
        maximum_balance_error=max(r['max_balance_error'] for r in audits),
        source_design_fallbacks=sum(r['source_fallback'] for r in audits),
        target_design_fallbacks=sum(r['target_fallback'] for r in audits),
        dispersion_failures_in_model=sum(r['dispersion_failure'] for r in exact),
        source_variance_ratio_range=[min(r['variance_ratio'] for r in exact), max(r['variance_ratio'] for r in exact)],
        note='Two declared CATE classes are distinct guarantees; nonlinear panel is outside the linear model.')
    (report/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    table = ['$p$ & $K$ & New length & Protected target & Normal target & Coverage & Paired change (MCSE)\\\\\\midrule']
    chosen = [r for r in rows if r['method']=='adaptive_universal' and r['pattern']=='weak'
              and ((r['profile']=='main' and r['source_count']==256) or (r['profile']=='boundary' and r['dimension']>=50))]
    for row in chosen:
        target = lookup[(row['profile'], row['cell_id'], 'target_universal')]
        normal = lookup[(row['profile'], row['cell_id'], 'target_normal')]
        change = next(r for r in paired if r['profile']==row['profile'] and r['cell_id']==row['cell_id'] and r['method']==row['method'])
        table.append('{} & {} & {:.4f} & {:.4f} & {:.4f} & {:.3f} & {:.4f} ({:.4f})\\\\'.format(
            row['dimension'], row['source_count'], row['length'], target['length'], normal['length'],
            row['coverage'], change['length_change'], change['length_change_mcse']))
    write_table(report/'continuous_table.tex', 'rrrrrrr', table)
    table = ['Pattern & New length & Target length & Coverage & Projection RMSE & Midpoint RMSE & Target RMSE\\\\\\midrule']
    for row in rows:
        if row['profile']=='boundary' and row['dimension']==20 and row['source_shift']==.4 and row['method']=='adaptive_universal':
            target = lookup[('boundary', row['cell_id'], 'target_universal')]
            table.append('{} & {:.4f} & {:.4f} & {:.3f} & {:.4f} & {:.4f} & {:.4f}\\\\'.format(
                row['pattern'].capitalize(), row['length'], target['length'], row['coverage'], row['rmse'],
                row['midpoint_rmse'], target['rmse']))
    write_table(report/'boundary_table.tex', 'lrrrrrr', table)
    table = ['$K$ & Naive, valid & Naive, weak & Protected weak & New length & Target length\\\\\\midrule']
    for count in [128, 512, 2048]:
        valid = next(r for r in audits if r['profile']=='stress' and r['source_count']==count and r['pattern']=='all_valid')
        weak = next(r for r in audits if r['profile']=='stress' and r['source_count']==count and r['pattern']=='weak')
        method = lookup[('stress', weak['cell_id'], 'adaptive_universal')]
        target = lookup[('stress', weak['cell_id'], 'target_universal')]
        table.append('{} & {:.3f} & {:.3f} & {:.3f} & {:.4f} & {:.4f}\\\\'.format(count,
            valid['naive_conditional_coverage'], weak['naive_conditional_coverage'], method['coverage'], method['length'], target['length']))
    write_table(report/'continuous_stress_table.tex', 'rrrrrr', table)
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
                         'axes.grid':True, 'grid.alpha':.2, 'pdf.fonttype':42})
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    for dimension, color in [(4, '#0072B2'), (10, '#D55E00'), (20, '#009E73')]:
        selection = sorted([r for r in rows if r['profile']=='main' and r['pattern']=='weak'
            and r['dimension']==dimension and r['method']=='adaptive_universal'], key=lambda r:r['source_count'])
        counts = [r['source_count'] for r in selection]
        for axis, baseline in zip(axes, ['target_universal', 'target_normal']):
            ratios = [r['length']/lookup[('main', r['cell_id'], baseline)]['length'] for r in selection]
            axis.plot(counts, ratios, 'o-', color=color, label='$p={}$'.format(dimension))
            axis.set_xscale('log', base=2); axis.set_xticks(counts); axis.set_xticklabels(counts)
            axis.set_xlabel('Number of sources K'); axis.axhline(1, color='black', ls=':', lw=1)
    axes[0].set_ylabel('Length / protected target-only'); axes[1].set_ylabel('Length / normal target-only')
    axes[0].legend(frameon=False); fig.tight_layout(); fig.savefig(str(report/'continuous_comparison.pdf')); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    selected = sorted([r for r in audits if r['profile']=='stress' and r['pattern']=='weak'], key=lambda r:r['source_count'])
    counts = [r['source_count'] for r in selected]
    axes[0].plot(counts, [r['naive_conditional_coverage'] for r in selected], 'o-', color='#D55E00', label='Naive conditional pooling')
    axes[0].plot(counts, [lookup[('stress', r['cell_id'], 'adaptive_universal')]['coverage'] for r in selected], 'o-', color='#0072B2', label='Protected population interval')
    axes[0].axhline(.95, color='black', ls=':'); axes[0].set_ylim(0,1.04); axes[0].set_ylabel('Coverage')
    for method, label, color in [('adaptive_universal','Universal range 2','#0072B2'),
                                 ('adaptive_declared_half','Declared range 0.5','#009E73'),
                                 ('target_universal','Protected target, range 2','#777777')]:
        axes[1].plot(counts, [lookup[('stress', r['cell_id'], method)]['length'] for r in selected], 'o-', label=label, color=color)
    axes[1].set_ylabel('Mean population interval length')
    for axis in axes:
        axis.set_xscale('log',base=2); axis.set_xticks(counts); axis.set_xticklabels(counts)
        axis.set_xlabel('Number of sources K'); axis.legend(fontsize=8,frameon=False)
    fig.tight_layout(); fig.savefig(str(report/'continuous_stress.pdf')); plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

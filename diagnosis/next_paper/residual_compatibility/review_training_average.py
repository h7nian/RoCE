"""Summarize paired reuse, independent confirmation, and stronger references."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_rows(path):
    return list(csv.DictReader(path.open()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    root = Path(parser.parse_args().root)
    report = root / 'report'
    report.mkdir(exist_ok=True)
    profiles = {name: read_rows(root/name/'metrics.csv') for name in ['paired','reused','confirmation','large_k']}
    target = read_rows(root/'target_baseline_audit/metrics.csv')
    target_lookup = {(r['profile'],r['cell_id'],r['method']):r for r in target}
    summary_rows = []
    for profile in ['confirmation','large_k']:
        for row in profiles[profile]:
            if row['method'] not in ['conditional_pilot','conditional_training_p004']:
                continue
            target_row = target_lookup[(profile,row['cell_id'],row['method'])]
            reference_length = float(row['mean_length']) / float(row['length_ratio_certified_target'])
            summary_rows.append(dict(profile=profile,source_count=row['source_count'],pattern=row['pattern'],
                repeats=row['repeats'],method=row['method'],coverage=row['coverage'],mean_length=row['mean_length'],
                rmse=row['rmse_midpoint'],bias=row['bias_midpoint'],target_aipw_length=reference_length,
                target_stratified_length=target_row['target_length'],target_stratified_rmse=target_row['target_rmse'],
                ratio_stratified=target_row['length_ratio']))
    with (report/'selected_metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summary_rows[0]));writer.writeheader();writer.writerows(summary_rows)
    table=['\\begin{tabular}{rllrrrr}\\toprule',
           '$K$ & Bias & Method & Coverage & Length & Stratified target & RMSE\\\\\\midrule']
    for row in summary_rows:
        if row['pattern']=='all_valid':continue
        method='Design' if row['method']=='conditional_pilot' else 'MGF (.04)'
        pattern='Weak' if row['pattern']=='weak_boundary' else 'Strong'
        table.append('{} & {} & {} & {:.3f} & {:.4f} & {:.4f} & {:.4f}\\\\'.format(
            row['source_count'],pattern,method,float(row['coverage']),float(row['mean_length']),
            float(row['target_stratified_length']),float(row['rmse'])))
    table.extend(['\\bottomrule','\\end{tabular}'])
    (report/'confirmation_table.tex').write_text('\n'.join(table)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.grid':True,'grid.alpha':.2,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    curves=[('paired','previous_corner','Previous corners, nE=250','#777777'),
            ('paired','conditional_pilot','Conditional design, nE=250','#0072B2'),
            ('reused','conditional_pilot','Conditional design, nE=500','#D55E00'),
            ('reused','conditional_training_p004','Training MGF (.04), nE=500','#009E73')]
    for profile,method,label,color in curves:
        rows=sorted([r for r in profiles[profile] if r['pattern']=='weak_boundary' and r['shared_scale']=='0.3' and r['method']==method],key=lambda r:int(r['source_count']))
        axes[0].plot([int(r['source_count']) for r in rows],[float(r['mean_length']) for r in rows],'o-',color=color,label=label)
    reference_row=next(r for r in profiles['reused'] if r['cell_id']=='16' and r['method']=='conditional_pilot')
    aipw=float(reference_row['mean_length'])/float(reference_row['length_ratio_certified_target'])
    stratified=float(target_lookup[('reused','16','conditional_pilot')]['target_length'])
    axes[0].axhline(aipw,color='black',ls='--',lw=1,label='Two-fold target AIPW')
    axes[0].axhline(stratified,color='black',ls=':',lw=1,label='Unsplit stratified target')
    axes[0].set_xscale('log',base=2);axes[0].set_xticks([8,32,128,512]);axes[0].set_xticklabels(['8','32','128','512'])
    axes[0].set_xlabel('Number of sources K');axes[0].set_ylabel('Mean interval length')
    axes[0].legend(frameon=True,facecolor='white',edgecolor='none',framealpha=.95,fontsize=7)
    for pattern,label,color in [('weak_boundary','Weak bias','#0072B2'),('strong','Strong bias','#D55E00')]:
        rows=sorted([r for r in profiles['confirmation'] if r['pattern']==pattern and r['method']=='conditional_pilot'],key=lambda r:int(r['source_count']))
        axes[1].plot([int(r['source_count']) for r in rows],[float(r['rmse_midpoint']) for r in rows],'o-',label=label,color=color)
    ref=target_lookup[('confirmation','3','conditional_pilot')]
    axes[1].axhline(float(ref['target_rmse']),color='black',ls=':',label='Stratified target RMSE')
    axes[1].set_xticks([128,512]);axes[1].set_xlabel('Number of sources K');axes[1].set_ylabel('Midpoint RMSE')
    axes[1].legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(str(report/'training_comparison.pdf'));plt.close(fig)
    differences=[]
    for cell in range(1,17):
        first=np.load(str(root/'paired/cells'/('cell_%03d.npz'%cell)))
        second=np.load(str(root/'reused/cells'/('cell_%03d.npz'%cell)))
        np.testing.assert_array_equal(first['target_reference_point'],second['target_reference_point'])
        for method in ['previous_corner','conditional_pilot','conditional_training_p004']:
            index=list(first['methods']).index(method)
            change=np.diff(second['interval'][:,index],axis=1)[:,0]-np.diff(first['interval'][:,index],axis=1)[:,0]
            differences.append(dict(cell_id=cell,method=method,mean_length_change=float(change.mean()),
                mcse=float(change.std(ddof=1)/np.sqrt(len(change)))))
    with (report/'paired_reuse_differences.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(differences[0]));writer.writeheader();writer.writerows(differences)
    failures=[]
    for profile in profiles:
        for path in (root/profile/'diagnostics').glob('*.json'):
            failures.extend(row['certificate_failure'] for row in json.loads(path.read_text()).values())
    summary=dict(main_settings=16,main_repeats=300,paired_reuse_same_data=True,
        confirmation_settings=4,confirmation_repeats=1000,exploratory_large_k_settings=3,large_k_repeats=300,
        total_distinct_setting_repetitions=9700,method_setting_rows=sum(len(v) for v in profiles.values()),
        observed_coverage_range=[min(float(r['coverage']) for rows in profiles.values() for r in rows),max(float(r['coverage']) for rows in profiles.values() for r in rows)],
        maximum_observed_certificate_failure=max(failures),
        actual_minimum_source_joint_probability=float(.35/(2*(1+np.exp(.6)))))
    (report/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    main()

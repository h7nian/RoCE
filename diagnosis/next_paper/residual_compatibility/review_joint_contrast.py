"""Review joint-score ablations and the exact-calibration diagnostic separately."""
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


def write_table(path,header,rows):
    lines=['\\begin{tabular}{'+header+'}\\toprule']+rows+['\\bottomrule','\\end{tabular}']
    path.write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('root');args=parser.parse_args()
    root=Path(args.root);report=root/'report';report.mkdir(exist_ok=True)
    primary={name:read_rows(root/name/'metrics.csv') for name in ['main','confirmation']}
    follow=read_rows(root/'followup/metrics.csv')
    balance=read_rows(root/'conditional_balance_fresh/metrics.csv')
    points=read_rows(root/'conditional_point_fresh/metrics.csv')
    stress=read_rows(root/'conditional_weak_stress/metrics.csv')
    follow_lookup={(r['profile'],r['cell_id'],r['method']):r for r in follow}
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.grid':True,'grid.alpha':.2,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    for method,label,color in [('separate_noise_separate_bias','Separate noise and drift','#777777'),
                              ('joint_noise_separate_bias','Joint noise only','#0072B2'),
                              ('separate_noise_joint_bias','Joint drift only','#D55E00'),
                              ('joint_noise_joint_bias','Joint noise and drift','#009E73')]:
        rows=sorted([r for r in primary['confirmation'] if r['pattern']=='weak_boundary' and r['method']==method],key=lambda r:int(r['source_count']))
        x=[int(r['source_count']) for r in rows]
        axes[0].plot(x,[float(r['source_length']) for r in rows],'o-',label=label,color=color)
        axes[1].plot(x,[float(r['mean_length']) for r in rows],'o-',label=label,color=color)
    target=float(follow_lookup[('confirmation','5','contrast_dispersion_design')]['target_range_length'])
    axes[1].axhline(target,color='black',ls=':',label='Range-adaptive target-only')
    for axis in axes:
        axis.set_xscale('log',base=2);axis.set_xticks([128,512,2048]);axis.set_xticklabels(['128','512','2048'])
        axis.set_xlabel('Number of sources K')
    axes[0].set_ylabel('Source-band length (before intersection)')
    axes[1].set_ylabel('Reported interval length')
    axes[0].legend(fontsize=7,frameon=True,edgecolor='none',facecolor='white')
    axes[1].legend(fontsize=7,frameon=True,edgecolor='none',facecolor='white')
    fig.tight_layout();fig.savefig(str(report/'joint_ablation.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    for method,label,color in [('heldout_adaptive','500 source outcomes','#0072B2'),('all_adaptive','1000 source outcomes','#009E73')]:
        rows=sorted([r for r in balance if r['pattern']=='weak_boundary' and r['method']==method],key=lambda r:int(r['source_count']))
        axes[0].plot([int(r['source_count']) for r in rows],[float(r['mean_length']) for r in rows],'o-',label=label,color=color)
    target=float(next(r for r in balance if r['method']=='all_adaptive')['target_length'])
    axes[0].axhline(target,color='black',ls=':',label='Unsmoothed target-only')
    axes[0].set_xscale('log',base=2);axes[0].set_xticks([128,512,2048]);axes[0].set_xticklabels(['128','512','2048'])
    axes[0].set_xlabel('Number of sources K');axes[0].set_ylabel('Mean interval length')
    axes[0].legend(fontsize=8,frameon=False)
    patterns=['weak_boundary','strong','arm_cancel'];locations=np.arange(3)
    chosen=[next(r for r in points if r['method']=='all_adaptive' and r['pattern']==p) for p in patterns]
    for offset,field,label,color in [(-.25,'target_rmse','Target-only','#777777'),(0,'rmse','Midpoint','#0072B2'),
                                    (.25,'projection_rmse','Target projection','#D55E00')]:
        axes[1].bar(locations+offset,[float(r[field]) for r in chosen],width=.24,label=label,color=color)
    axes[1].set_xticks(locations);axes[1].set_xticklabels(['Weak bias','Strong bias','Arm cancellation'])
    axes[1].set_ylabel('Population TATE RMSE');axes[1].legend(fontsize=8,frameon=False)
    fig.tight_layout();fig.savefig(str(report/'conditional_balance.pdf'));plt.close(fig)
    rows=['$K$ & Source pattern & 500 outcomes & 1000 outcomes & Target-only & Coverage\\\\\\midrule']
    for cell in range(1,6):
        pair=[next(r for r in balance if r['cell_id']==str(cell) and r['method']==m) for m in ['heldout_adaptive','all_adaptive']]
        label={'weak_boundary':'Weak','strong':'Strong','arm_cancel':'Cancellation'}[pair[0]['pattern']]
        rows.append('{} & {} & {:.4f} & {:.4f} & {:.4f} & {:.3f}\\\\'.format(pair[0]['source_count'],label,
            float(pair[0]['mean_length']),float(pair[1]['mean_length']),float(pair[1]['target_length']),float(pair[1]['coverage'])))
    write_table(report/'balance_table.tex','rlrrrr',rows)
    rows=['Source pattern & Midpoint RMSE & Projection RMSE & Target RMSE & Projection bias\\\\\\midrule']
    for r in chosen:
        label={'weak_boundary':'Weak','strong':'Strong','arm_cancel':'Cancellation'}[r['pattern']]
        rows.append('{} & {:.4f} & {:.4f} & {:.4f} & {:.4f}\\\\'.format(label,float(r['rmse']),float(r['projection_rmse']),float(r['target_rmse']),float(r['projection_bias'])))
    write_table(report/'point_table.tex','lrrrr',rows)
    rows=['$K$ & All-valid naive coverage & Weak naive coverage & Protected weak coverage & Protected length\\\\\\midrule']
    variance_audit=[]
    for count in [128,512,2048]:
        valid=next(r for r in stress if int(r['source_count'])==count and r['pattern']=='all_valid' and r['method']=='all_adaptive')
        weak=next(r for r in stress if int(r['source_count'])==count and r['pattern']=='weak_boundary' and r['method']=='all_adaptive')
        rows.append('{} & {:.3f} & {:.3f} & {:.3f} & {:.4f}\\\\'.format(count,float(valid['naive_coverage']),
            float(weak['naive_coverage']),float(weak['coverage']),float(weak['mean_length'])))
        packet=np.load(str(root/'conditional_weak_stress/weak_stress'/('cell_%03d.npz'%int(valid['cell_id']))))
        index=list(packet['methods']).index('all_adaptive')
        interval=packet['naive_interval'][:,index]
        if np.any(abs(interval)>=1):raise AssertionError('Cannot reconstruct a clipped normal radius')
        estimated_variance=((interval[:,1]-interval[:,0])/(2*1.959963984540054))**2
        observed_variance=np.var(packet['source_center'][:,index],ddof=1)
        variance_audit.append(dict(source_count=count,ratio_mean_estimated_to_empirical_variance=float(estimated_variance.mean()/observed_variance)))
    write_table(report/'weak_stress_table.tex','rrrrr',rows)
    (report/'naive_variance_audit.json').write_text(json.dumps(variance_audit,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(9,3.6))
    weak=sorted([r for r in stress if r['pattern']=='weak_boundary' and r['method']=='all_adaptive'],key=lambda r:int(r['source_count']))
    x=[int(r['source_count']) for r in weak]
    axes[0].plot(x,[float(r['naive_coverage']) for r in weak],'o-',color='#D55E00',label='Unprotected pooling')
    axes[0].plot(x,[float(r['coverage']) for r in weak],'o-',color='#0072B2',label='Conditional compatibility')
    axes[0].axhline(.95,color='black',ls=':',lw=1);axes[0].set_ylim(0,1.04);axes[0].set_ylabel('Coverage')
    axes[0].legend(fontsize=8,frameon=False)
    axes[1].plot(x,[float(r['mean_length']) for r in weak],'o-',color='#0072B2',label='Conditional compatibility')
    axes[1].plot(x,[float(r['target_length']) for r in weak],'o:',color='black',label='Protected target-only')
    axes[1].set_ylabel('Mean interval length');axes[1].legend(fontsize=8,frameon=False)
    for axis in axes:
        axis.set_xscale('log',base=2);axis.set_xticks(x);axis.set_xticklabels([str(v) for v in x]);axis.set_xlabel('Number of sources K')
    fig.tight_layout();fig.savefig(str(report/'weak_bias_stress.pdf'));plt.close(fig)
    differences=[]
    for profile in primary:
        cells=sorted(set(int(r['cell_id']) for r in primary[profile]))
        for cell in cells:
            data=np.load(str(root/profile/'cells'/('cell_%03d.npz'%cell)))
            names=list(data['methods']);base=names.index('separate_noise_separate_bias')
            for method in ['joint_noise_separate_bias','separate_noise_joint_bias','joint_noise_joint_bias']:
                index=names.index(method)
                change=np.diff(data['interval'][:,index],axis=1)[:,0]-np.diff(data['interval'][:,base],axis=1)[:,0]
                differences.append(dict(profile=profile,cell_id=cell,method=method,mean_change=float(change.mean()),
                    paired_mcse=float(change.std(ddof=1)/np.sqrt(len(change)))))
    with (report/'paired_joint_differences.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(differences[0]));writer.writeheader();writer.writerows(differences)
    all_rows=primary['main']+primary['confirmation']+follow+balance+points+stress+read_rows(root/'conditional_balance_paired/metrics.csv')
    summary=dict(primary_main_repetitions=4500,primary_confirmation_repetitions=6000,
        exact_balance_fresh_repetitions=2500,point_rule_fresh_repetitions=3000,weak_stress_fresh_repetitions=3000,
        distinct_reported_setting_repetitions=19000,paired_reanalyses_not_recounted=True,
        method_setting_rows=len(all_rows),observed_coverage_range=[min(float(r['coverage']) for r in all_rows),max(float(r['coverage']) for r in all_rows)],
        note='Exact-balance score changes and joint-score inference are distinct procedures. No high-dimensional RoCE theorem is claimed.')
    (report/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()

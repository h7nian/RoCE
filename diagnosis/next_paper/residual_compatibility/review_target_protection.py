"""Audit and report paired target protection, model error, and design stability."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def table(path,columns,rows):
    path.write_text('\n'.join(['\\begin{tabular}{'+columns+'}\\toprule']+rows+['\\bottomrule','\\end{tabular}'])+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path);args=parser.parse_args()
    root=args.root;report=root/'report';report.mkdir(exist_ok=True)
    required=['target','misspecification_v2','misspecification_bounded','projection_stress','projection_stress_bounded','design_stability']
    if any(not (root/p/'COMPLETE').exists() for p in required):raise RuntimeError('Finish all prescribed studies before reporting')
    rows=[];paired=[];audits=[];parity=[]
    for profile in ['target','misspecification_bounded']:
        for path in sorted((root/profile).glob('cell_*.json')):
            packet=json.loads(path.read_text());data=np.load(str(path.with_suffix('.npz')));names=packet['methods']
            if profile=='misspecification_bounded':
                old_path=root/'misspecification_v2'/path.name;old=json.loads(old_path.read_text());old_data=np.load(str(old_path.with_suffix('.npz')))
                for i,method in enumerate(old['methods']):
                    diff=float(np.max(np.abs(data['interval'][:,names.index(method)]-old_data['interval'][:,i])))
                    if diff>1e-12:raise AssertionError('Paired extension changed a prior interval')
                parity.append(dict(profile=profile,cell_id=packet['setting']['cell_id'],max_difference=0.))
            for i,metric in enumerate(packet['metrics']):
                interval=data['interval'][:,i];conditional=data['conditional'][:,i];truth=data['conditional_truth']
                covered=(interval[:,0]<=.1)&(.1<=interval[:,1]);cond=(conditional[:,0]<=truth)&(truth<=conditional[:,1])
                length=np.diff(interval,axis=1)[:,0]
                for key,value in dict(coverage=covered.mean(),conditional_coverage=cond.mean(),length=length.mean()).items():
                    if abs(metric[key]-value)>1e-12:raise AssertionError('Metric audit failed')
                count=int(covered.sum());total=len(covered)
                metric.update(profile=profile,length_mcse=float(length.std(ddof=1)/np.sqrt(total)),
                    coverage_lower=float(beta.ppf(.025,count,total-count+1)) if count else 0.,
                    coverage_upper=float(beta.ppf(.975,count+1,total-count)) if count<total else 1.)
                rows.append(metric)
            source_index=names.index('source_universal')
            for method in ['source_range','source_bernstein','source_adaptive']:
                i=names.index(method)
                np.testing.assert_allclose(data['conditional'][:,i],data['conditional'][:,source_index],rtol=0,atol=1e-14)
                change=np.diff(data['interval'][:,i],axis=1)[:,0]-np.diff(data['interval'][:,source_index],axis=1)[:,0]
                paired.append(dict(profile=profile,cell_id=packet['setting']['cell_id'],method=method,
                    mean_length_change=float(change.mean()),paired_mcse=float(change.std(ddof=1)/np.sqrt(len(change)))))
            diagnostics=packet['diagnostics']
            audits.append(dict(profile=profile,**packet['setting'],repeats=len(diagnostics),
                certificate_failures=sum(r['certificate_failed'] for r in diagnostics),
                range_mean=float(np.mean([r['range_upper'] for r in diagnostics])),
                range_cap_fraction=float(np.mean([r['range_upper']>=2-1e-12 for r in diagnostics])),
                variance_upper_mean=float(np.mean([r['variance_upper'] for r in diagnostics])),
                source_failures=sum(r['source_failed'] for r in diagnostics),target_failures=sum(r['target_failed'] for r in diagnostics)))
    projection=[]
    for path in sorted((root/'projection_stress_bounded').glob('k*.json')):
        packet=json.loads(path.read_text());data=np.load(str(path.with_suffix('.npz')))
        old=np.load(str(root/'projection_stress'/path.with_suffix('.npz').name))
        diff=float(np.max(np.abs(data['conditional'][:,:2]-old['conditional'])))
        if diff>1e-12:raise AssertionError('Fixed-design paired extension changed its controls')
        parity.append(dict(profile='projection_stress',cell_id=path.stem,max_difference=diff))
        for i,row in enumerate(packet['metrics']):
            conditional=data['conditional'][:,i];pop=data['population'][:,i]
            if abs(row['conditional_coverage']-np.mean((conditional[:,0]<=.1)&(.1<=conditional[:,1])))>1e-12:raise AssertionError('Conditional coverage mismatch')
            if abs(row['population_length']-np.diff(pop,axis=1).mean())>1e-12:raise AssertionError('Population width mismatch')
            projection.append(row)
    designs=[json.loads(p.read_text()) for p in sorted((root/'design_stability').glob('cell_*.json'))]
    for packet in designs:
        observed=np.mean([r['any_failure'] for r in packet['records']])
        if abs(observed-packet['metrics']['any_source_failure'])>1e-12:raise AssertionError('Design failure audit failed')
    for filename,values in [('metrics.csv',rows),('paired_target_changes.csv',paired),('certificate_audit.csv',audits),
                            ('projection_metrics.csv',projection),('design_metrics.csv',[p['metrics'] for p in designs])]:
        fields=sorted({key for row in values for key in row})
        with (report/filename).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(values)
    (report/'paired_parity.json').write_text(json.dumps(parity,indent=2)+'\n')
    lookup={(r['profile'],r['cell_id'],r['method']):r for r in rows}
    selected=[r for r in rows if r['profile']=='target' and r['method']=='source_bernstein'
              and r['pattern']=='weak' and r['source_count']>=256]
    tab=['$p$ & $K$ & $h$ & Universal & Range & Bernstein & Target EB & Target normal & Coverage\\\\\\midrule']
    for r in selected:
        related=[lookup[('target',r['cell_id'],m)] for m in ['source_universal','source_range','source_bernstein','target_bernstein','target_normal']]
        tab.append('{} & {} & {:.2f} & {} & {:.3f}\\\\'.format(r['dimension'],r['source_count'],r['heterogeneity'],
            ' & '.join('{:.4f}'.format(q['length']) for q in related),r['coverage']))
    table(report/'target_table.tex','rrrrrrrrr',tab)
    tab=['$\\eta$ & Old conditional & Old population & Envelope conditional & Bounded conditional & Bounded width\\\\\\midrule']
    for eta in [0.,.1,.2,.4]:
        q=[next(r for r in projection if r['source_count']==262144 and r['nonlinearity']==eta and r['method']==m)
           for m in ['uncorrected','projection_repaired','bounded_invalid']]
        tab.append('{:.1f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.5f}\\\\'.format(eta,
            q[0]['conditional_coverage'],q[0]['population_coverage'],q[1]['conditional_coverage'],
            q[2]['conditional_coverage'],q[2]['conditional_length']))
    table(report/'projection_table.tex','rrrrrr',tab)
    tab=['$n_s$ & $p$ & $K$ & Strength & Any failure & Eligible after removal & Sufficient bound\\\\\\midrule']
    for packet in designs:
        r=packet['metrics']
        if r['cell_id'] in [8,10,12,24,36,37,38,39]:
            tab.append('{} & {} & {} & {:.1f} & {:.3f} & {:.3f} & {:.3g}\\\\'.format(
                r['source_size'],r['dimension'],r['source_count'],r['overlap_strength'],r['any_source_failure'],
                r['eligible_after_count_adjustment'],r['any_source_gram_failure']))
    table(report/'design_table.tex','rrrrrrr',tab)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.grid':True,'grid.alpha':.2,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(9,3.6))
    for axis,dimension in zip(axes,[4,10]):
        selected=sorted([r for r in rows if r['profile']=='target' and r['dimension']==dimension and r['heterogeneity']==0
                         and r['method']=='source_bernstein'],key=lambda r:r['source_count'])
        x=[r['source_count'] for r in selected]
        for method,label,color,style in [('source_universal','Universal range','#D55E00','--'),
            ('source_bernstein','Estimated variance / range','#0072B2','-'),
            ('target_bernstein','Improved protected target','#777777','--'),('target_normal','Normal target','#222222',':')]:
            axis.plot(x,[lookup[('target',r['cell_id'],method)]['length'] for r in selected],marker='o',ls=style,color=color,label=label)
        axis.set_xscale('log',base=2);axis.set_xticks(x);axis.set_xticklabels(x);axis.set_xlabel('Number of sources K')
        axis.set_ylabel('Mean population interval length, p={}'.format(dimension))
    axes[0].legend(fontsize=8,frameon=False);fig.tight_layout();fig.savefig(str(report/'target_protection.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.6))
    for method,label,color in [('uncorrected','Uncorrected','#D55E00'),('projection_repaired','Declared envelope','#0072B2'),('bounded_invalid','Bounded invalid means','#009E73')]:
        selected=sorted([r for r in projection if r['nonlinearity']==.4 and r['method']==method],key=lambda r:r['source_count'])
        x=[r['source_count'] for r in selected]
        axes[0].plot(x,[r['conditional_coverage'] for r in selected],'o-',color=color,label=label)
        axes[1].plot(x,[r['conditional_length'] for r in selected],'o-',color=color,label=label)
    axes[0].axhline(.95,color='black',ls=':');axes[0].set_ylim(-.03,1.04);axes[0].set_ylabel('Conditional coverage')
    axes[1].set_yscale('log');axes[1].set_ylabel('Conditional interval length')
    for axis in axes:
        axis.set_xscale('log',base=2);axis.set_xticks(x);axis.set_xticklabels(x,rotation=20);axis.set_xlabel('Number of sources K')
    axes[0].legend(fontsize=8,frameon=False);fig.tight_layout();fig.savefig(str(report/'projection_stress.pdf'));plt.close(fig)
    summary=dict(outcome_setting_repetitions=7200,design_only_repetitions=7800,
        paired_extensions_not_recounted=4800,reported_outcome_settings=36,design_settings=len(designs),
        target_certificate_failures=sum(r['certificate_failures'] for r in audits),
        source_design_failures=sum(r['source_failures'] for r in audits),
        target_design_failures=sum(r['target_failures'] for r in audits),
        paired_parity_checks=len(parity),maximum_paired_difference=max(r['max_difference'] for r in parity),
        note='Fixed-design and random-design experiments are distinct; invalid eta=.2 random-DGP attempt is preserved separately.')
    (report/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()

"""Endpoint audits, paired uncertainty and publication assets for v22."""
import argparse
import json
from pathlib import Path

import numpy as np
from review_design_followup import write_csv, binomial_interval
from weak_bias_precision import local_length_lower_bound


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();root=args.root;report=root/'report';report.mkdir(exist_ok=True)
    for panel in ['budget','local_weak_bias']:
        if not (root/panel/'COMPLETE').exists():raise ValueError('Panel incomplete: '+panel)
    rows=[];pairs=[]
    for path in sorted((root/'budget').glob('*_*.json')):
        d=json.loads(path.read_text())
        if 'metrics' not in d:continue
        data=np.load(str(path.with_suffix('.npz')));names=d['methods'];n=len(data['interval'])
        lengths=np.diff(data['interval'],axis=-1)[...,0]
        for j,name in enumerate(names):
            interval=data['interval'][:,j];conditional=data['conditional'][:,j]
            covered=(interval[:,0]<=.1)&(.1<=interval[:,1])
            conditional_coverage=np.mean((conditional[:,0]<=data['conditional_truth'])&(data['conditional_truth']<=conditional[:,1]))
            old=d['metrics'][j]
            np.testing.assert_allclose([old['length'],old['coverage'],old['conditional_coverage']],
                [lengths[:,j].mean(),covered.mean(),conditional_coverage],rtol=0,atol=1e-12)
            lo,hi=binomial_interval(int(covered.sum()),n)
            rows.append(dict(old,length_mcse=float(lengths[:,j].std(ddof=1)/np.sqrt(n)),coverage_lower=lo,coverage_upper=hi))
        for method in ['conditional_reallocation','shared_composition','guard2','guard4','guard8','precision_gate','guard4_gate']:
            for comparator in ['legacy','target_protected','target_normal']:
                j,k=names.index(method),names.index(comparator);diff=lengths[:,j]-lengths[:,k]
                pairs.append(dict(d['setting'],method=method,comparator=comparator,repeats=n,
                    difference=float(diff.mean()),paired_mcse=float(diff.std(ddof=1)/np.sqrt(n)),
                    length_ratio=float(lengths[:,j].mean()/lengths[:,k].mean())))
    write_csv(report/'budget_metrics.csv',rows);write_csv(report/'budget_pairs.csv',pairs)
    local=[]
    for path in sorted((root/'local_weak_bias').glob('k*.json')):
        d=json.loads(path.read_text());data=np.load(str(path.with_suffix('.npz')))
        for j,row in enumerate(d['metrics']):
            interval=data['interval'][:,j];covered=(interval[:,0]<=data['truth'])&(data['truth']<=interval[:,1])
            lengths=np.diff(interval,axis=1)[:,0]
            np.testing.assert_allclose([row['coverage'],row['length']],[covered.mean(),lengths.mean()],rtol=0,atol=1e-12)
            lo,hi=binomial_interval(int(covered.sum()),len(covered))
            local.append(dict(row,length_mcse=float(lengths.std(ddof=1)/np.sqrt(len(lengths))),coverage_lower=lo,coverage_upper=hi))
    write_csv(report/'local_metrics.csv',local)
    bounds=[local_length_lower_bound(n,nt,k,c,experiment=mode) for mode in ['fixed_arm','iid_randomized']
            for n in [125,500,1000,2000] for nt in [n,4*n] for k in [16,64,256,1024,4096,16384] for c in [.1,.2]]
    (report/'local_lower_bounds.json').write_text(json.dumps(bounds,indent=2)+'\n')
    def cell(panel,index):return {r['method']:r for r in rows if r['panel']==panel and r['cell_id']==index}
    lines=[]
    for index in range(1,5):
        c=cell('confirmation',index);r=c['shared_composition']
        pair=next(x for x in pairs if x['panel']=='confirmation' and x['cell_id']==index
                  and x['method']=='shared_composition' and x['comparator']=='legacy')
        lines.append('%d & %d & %.2f & %s & %.4f & %.4f & %.4f & %.3f & %.4f (%.6f) \\\\' %
            (r['dimension'],r['source_count'],r['heterogeneity'],r['pattern'].replace('_',' '),
             c['legacy']['length'],r['length'],c['target_normal']['length'],r['coverage'],pair['difference'],pair['paired_mcse']))
    (report/'budget_table.tex').write_text('\\begin{tabular}{rrrlrrrrr}\n\\toprule\n'
        '$p$ & $K$ & $h$ & Pattern & Legacy & Shared budget & Target Wald & Coverage & Change (MCSE) \\\\\n\\midrule\n'
        +'\n'.join(lines)+'\n\\bottomrule\n\\end{tabular}\n')
    lines=[]
    for index in [2,4,6,8,9,10]:
        c=cell('removal',index);r=c['precision_gate']
        lines.append('%d & %d & %d & %.3f & %.4f & %.4f & %.4f & %.4f & %.3f \\\\' %
            (r['source_size'],r['dimension'],r['source_count'],r['singular_fraction'],
             c['legacy']['length'],c['shared_composition']['length'],c['guard4']['length'],r['length'],r['fallback_rate']))
    (report/'design_budget_table.tex').write_text('\\begin{tabular}{rrrrrrrrr}\n\\toprule\n'
        '$n_s$ & $p$ & $K$ & Singular fraction & Legacy & Shared & Guard 4 & Gate & Gate abstention \\\\\n\\midrule\n'
        +'\n'.join(lines)+'\n\\bottomrule\n\\end{tabular}\n')
    lines=[]
    for count in [16,256,4096,16384]:
        chosen={r['method']:r for r in local if r['source_count']==count and r['pattern']=='moment_matched' and r['shift_constant']==.25}
        s,t,naive=[chosen[k] for k in ['source_protected','target_wald','naive_source_wald']]
        lines.append('%d & %.4f & %.3f & %.4f & %.3f & %.4f & %.3f \\\\' %
            (count,s['length'],s['coverage'],t['length'],t['coverage'],naive['length'],naive['coverage']))
    (report/'local_table.tex').write_text('\\begin{tabular}{rrrrrrr}\n\\toprule\n'
        '$K$ & Protected length & Coverage & Target Wald length & Coverage & Naive length & Coverage \\\\\n\\midrule\n'
        +'\n'.join(lines)+'\n\\bottomrule\n\\end{tabular}\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.grid':True,'grid.alpha':.2,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for method,label,color in [('source_protected','Protected source-assisted','#0072B2'),
            ('target_wald','Target-only (Wald)','#666666'),('naive_source_wald','Naive source pooling','#D55E00')]:
        chosen=sorted([r for r in local if r['method']==method and r['pattern']=='moment_matched'
                       and r['shift_constant']==.25],key=lambda x:x['source_count'])
        counts=[r['source_count'] for r in chosen]
        axes[0].semilogx(counts,[r['coverage'] for r in chosen],'-o',color=color,label=label)
        axes[1].loglog(counts,[r['length'] for r in chosen],'-o',color=color,label=label)
    axes[0].axhline(.95,color='black',ls=':',lw=1)
    axes[0].set_ylim(0,1.03);axes[0].set_ylabel('Observed coverage')
    axes[1].set_ylabel('Mean confidence interval length')
    for ax in axes:ax.set_xlabel('Number of sources K')
    axes[0].legend(frameon=False,fontsize=8);fig.tight_layout()
    fig.savefig(str(report/'local_weak_bias.pdf'));fig.savefig(str(report/'local_weak_bias.png'),dpi=160);plt.close(fig)
    (report/'summary.json').write_text(json.dumps(dict(fresh_setting_repetitions=4800,
        paired_setting_reanalyses=2800,budget=rows,pairs=pairs,local=local),indent=2)+'\n')
    print('Verified',len(rows),'budget rows,',len(pairs),'paired comparisons and',len(local),'local rows')


if __name__=='__main__':main()

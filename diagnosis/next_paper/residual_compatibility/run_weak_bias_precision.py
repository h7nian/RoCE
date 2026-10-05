"""Exact intercept-model experiments on the moment-matched weak-bias scale."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from continuous_balance import linear_balance_design, conditional_linear_bands
from weak_bias_precision import draw_valid_flags, local_length_lower_bound


def intercept_summary(successes, arm_size, design):
    values=np.asarray(successes,dtype=float)
    return dict(estimate=values/arm_size,
        variance=values*(arm_size-values)/(arm_size**2*(arm_size-1)),
        **{key:np.full(values.shape,design[key]) for key in ['weight_square_sum','diagonal_square_sum',
            'correction_frobenius','max_diagonal','balance_error','negative_weight_fraction']})


def run_cell(source_count, pattern, constant, repeats, seed, arm_size=500):
    design=linear_balance_design(np.ones((arm_size,1)),np.ones(1))
    shift=constant/np.sqrt(arm_size)/source_count**.25
    labels=['source_protected','target_protected','target_wald','naive_source_wald']
    records=[];truths=[];valid_counts=[]
    for repeat in range(repeats):
        rng=np.random.RandomState(seed+repeat)
        flags=draw_valid_flags(rng,source_count)
        theta=.5+shift if pattern=='moment_matched' else .5
        probabilities=np.full(source_count,theta)
        if pattern=='moment_matched':probabilities[~flags]=.5-4*shift
        if pattern=='same_direction':probabilities[~flags]=.5+5*shift
        target_values=[intercept_summary(rng.binomial(arm_size,p),arm_size,design) for p in [theta,.4]]
        source_values=[intercept_summary(rng.binomial(arm_size,p),arm_size,design)
                       for p in [probabilities,np.full(source_count,.4)]]
        target={key:np.array([v[key] for v in target_values]) for key in target_values[0]}
        source={key:np.stack([v[key] for v in source_values],axis=-1) for key in source_values[0]}
        result=conditional_linear_bands(source,target,int(np.ceil(.75*source_count)),
            source_failure=.0125,target_failure=.0125,dispersion_failure=.025)
        target_point=target['estimate'][0]-target['estimate'][1]
        target_radius=np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/.05))
        source_point=float(np.mean(source['estimate'][:,0]-source['estimate'][:,1]))
        normal=norm.isf(.025)
        intervals=[result['conditional_interval'],target_point+target_radius*np.array([-1.,1.]),
            target_point+normal*np.sqrt(target['variance'].sum())*np.array([-1.,1.]),
            source_point+normal*np.sqrt(source['variance'].sum()/source_count**2)*np.array([-1.,1.])]
        records.append(intervals);truths.append(theta-.4)
        valid_counts.append(source_count if pattern=='all_valid' else int(flags.sum()))
    arrays=dict(interval=np.clip(np.array(records),-1,1),truth=np.array(truths),valid_count=np.array(valid_counts))
    metrics=[]
    for j,name in enumerate(labels):
        interval=arrays['interval'][:,j]
        metrics.append(dict(source_count=source_count,pattern=pattern,shift_constant=constant,
            source_arm_size=arm_size,method=name,repeats=repeats,
            coverage=float(np.mean((interval[:,0]<=arrays['truth'])&(arrays['truth']<=interval[:,1]))),
            length=float(np.diff(interval,axis=1).mean()),
            scaled_length=float(np.diff(interval,axis=1).mean()*np.sqrt(arm_size)*source_count**.25)))
    return arrays,metrics


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--repeats',type=int,default=200)
    args=parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats<2:
        raise ValueError('Use FACE-HD scratch and at least two repeats')
    output=args.root/'local_weak_bias';output.mkdir()
    counts=[16,64,256,1024,4096,16384]
    patterns=[('all_valid',0.),('moment_matched',.1),('moment_matched',.25),('same_direction',.25)]
    (output/'configuration.json').write_text(json.dumps(dict(source_counts=counts,patterns=patterns,
        repeats=args.repeats,n_per_site=1000,seed_offset=33190000,
        scope='Balanced fixed treatment groups; exact intercept Bernoulli model; no composition enlargement',
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}),indent=2)+'\n')
    bounds=[local_length_lower_bound(n,nt,k,c) for n in [125,500,2000]
            for nt in [n,4*n] for k in counts for c in [.1,.2]]
    (output/'exact_lower_bounds.json').write_text(json.dumps(bounds,indent=2)+'\n')
    rows=[]
    for index,count in enumerate(counts):
        for pattern,constant in patterns:
            arrays,metrics=run_cell(count,pattern,constant,args.repeats,33190000+10000*index)
            stem=output/('k%d_%s_c%03d'%(count,pattern,int(100*constant)))
            np.savez_compressed(str(stem)+'.npz',**arrays)
            stem.with_suffix('.json').write_text(json.dumps(dict(methods=[m['method'] for m in metrics],metrics=metrics),indent=2)+'\n')
            rows.extend(metrics);print('Local weak bias',count,pattern,constant,flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'COMPLETE').write_text('Every outcome repetition retained; latent source validity is used only for generating data.\n')


if __name__=='__main__':main()

"""Fixed-design Bernoulli stress isolating nonlinear leave-out variance bias."""
import argparse
import csv
import json
from pathlib import Path
import time
import numpy as np

from continuous_balance import linear_balance_design, conditional_linear_bands
from population_protection import projection_geometry, projection_error_budget, expand_conditional_interval


def fixed_geometry(arm_size=500):
    if arm_size % 4:
        raise ValueError('Arm size must be divisible by four')
    values=np.array([-1.,0.,1.]); counts=np.array([arm_size//4,arm_size//2,arm_size//4])
    x=np.repeat(values,counts); matrix=np.column_stack((np.ones(arm_size),x))
    design=linear_balance_design(matrix,np.array([1.,0.]))
    group_diagonal=np.array([design['diagonal'][np.flatnonzero(x==value)[0]] for value in values])
    constants={name:design[name] for name in ['weight_square_sum','diagonal_square_sum','correction_frobenius','max_diagonal','balance_error','negative_weight_fraction']}
    constants.update(projection_geometry(design))
    return values,counts,group_diagonal,constants


def summaries(successes,values,counts,diagonal,constants):
    total=counts.sum(); mean=successes.sum(axis=-1)/total
    slope=successes.dot(values)/np.sum(counts*values**2)
    fitted=mean[...,None]+slope[...,None]*values
    variance=np.sum(diagonal*(.5*successes-fitted*(successes-.5*counts)),axis=-1)
    return dict(estimate=mean,variance=variance,**{name:np.broadcast_to(value,mean.shape) for name,value in constants.items()})


def run_cell(count,eta,repeats,seed,bounded_invalid=False):
    values,counts,diagonal,constants=fixed_geometry()
    nonlinear=values**2-.5
    weak_shift=2*np.sqrt(.5/1000)*count**(-.25)
    probabilities=np.full((count,3),.5)
    probabilities[:count//4]+=weak_shift+eta*nonlinear
    template={name:np.full(count,value) for name,value in constants.items()}
    template['estimate']=np.zeros(count)
    budget=projection_error_budget(template,.2*np.sqrt(500),count//4)
    budgets={'mu1':budget,'mu0':dict(variance_bias_sum=0.,weighted_projection_error=0.),'contrast':budget}
    corrections=[None,budgets]
    labels=['uncorrected','projection_repaired']
    if bounded_invalid:
        universal=projection_error_budget(template,.5*np.sqrt(500),count//4)
        corrections.append({'mu1':universal,'mu0':universal,'contrast':{key:2*value for key,value in universal.items()}})
        labels.append('bounded_invalid')
    expectation_variance_bias=np.sum(counts*diagonal*(weak_shift+eta*nonlinear)*eta*nonlinear)
    true_dispersion=.25*.75*weak_shift**2
    mean_statistic=true_dispersion-(count-1)/count**2*(count//4)*expectation_variance_bias
    records=[]; started=time.time()
    for repeat in range(repeats):
        # Same stream within K across amplitudes, including identical target draws.
        rng=np.random.RandomState(seed+repeat)
        target_t=rng.binomial(counts,.5);target_c=rng.binomial(counts,.4)
        target_arms=[summaries(v,values,counts,diagonal,constants) for v in [target_t,target_c]]
        target={key:np.array([arm[key] for arm in target_arms]) for key in target_arms[0]}
        treated=rng.binomial(counts,probabilities)
        control=rng.binomial(counts,np.full((count,3),.4))
        arms=[summaries(v,values,counts,diagonal,constants) for v in [treated,control]]
        source={key:np.stack([arm[key] for arm in arms],axis=-1) for key in arms[0]}
        values_out=[]
        for correction in corrections:
            result=conditional_linear_bands(source,target,3*count//4,dispersion_budgets=correction)
            conditional=result['conditional_interval']
            population=expand_conditional_interval(conditional,2*np.sqrt(np.log(200.)/2000))
            values_out.append(dict(conditional=conditional,population=population,
                source_noise=result['source_noise'],dispersion=result['dispersion_radius'],
                source_center=result['source_center'],source_band=result['source_band']))
        records.append(values_out)
    arrays={name:np.array([[value[name] for value in r] for r in records])
            for name in ['conditional','population','source_noise','dispersion','source_center','source_band']}
    metrics=[]
    for j,name in enumerate(labels):
        cond=arrays['conditional'][:,j];pop=arrays['population'][:,j]
        metrics.append(dict(source_count=count,nonlinearity=eta,method=name,repeats=repeats,
            conditional_coverage=float(np.mean((cond[:,0]<=.1)&(.1<=cond[:,1]))),
            population_coverage=float(np.mean((pop[:,0]<=.1)&(.1<=pop[:,1]))),
            conditional_length=float(np.diff(cond,axis=1).mean()),population_length=float(np.diff(pop,axis=1).mean()),
            zero_dispersion_fraction=float(np.mean(arrays['dispersion'][:,j]==0))))
    audit=dict(source_count=count,nonlinearity=eta,weak_shift=weak_shift,
               true_dispersion=true_dispersion,expected_uncorrected_dispersion=mean_statistic,
               variance_centering_allowance=(count-1)/count**2*budget['variance_bias_sum'],
               runtime_seconds=time.time()-started)
    return arrays,metrics,audit


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--repeats',type=int,default=200);parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--bounded-invalid',action='store_true');parser.add_argument('--run-name')
    args=parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats<2:
        raise ValueError('Use FACE-HD scratch and at least two repetitions')
    name=args.run_name or ('projection_pilot' if args.pilot else 'projection_stress')
    if Path(name).name!=name or name in ['.','..']:raise ValueError('Use one output-directory name')
    output=args.root/name;output.mkdir()
    seed_offset=27190000 if args.pilot else 28190000
    counts=[128,2048,32768,262144] if not args.pilot else [128,262144]
    amplitudes=[0.,.1,.2,.4] if not args.pilot else [0.,.4]
    (output/'configuration.json').write_text(json.dumps(dict(source_counts=counts,amplitudes=amplitudes,
        repeats=args.repeats,n_per_site=1000,seed_offset=seed_offset,declared_treated_envelope=.2,bounded_invalid=args.bounded_invalid,
        scope='Conditional fixed-design counterexample; target CATE is constant, all valid means linear; population enlargement is retained for diagnosis.'),indent=2)+'\n')
    rows=[]
    for index,count in enumerate(counts):
        for eta in amplitudes:
            arrays,metrics,audit=run_cell(count,eta,args.repeats,seed_offset+index*10000,args.bounded_invalid)
            stem=output/('k{}_eta{}'.format(count,int(round(100*eta))))
            np.savez_compressed(str(stem)+'.npz',**arrays)
            Path(str(stem)+'.json').write_text(json.dumps(dict(metrics=metrics,audit=audit),indent=2)+'\n')
            rows.extend(metrics);print('Projection stress',count,eta,flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'COMPLETE').write_text('All conditional and population intervals retained.\n')


if __name__=='__main__':main()

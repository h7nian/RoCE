"""Prespecified paired TATE noise/drift ablations, preserving every repetition."""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
import numpy as np
from fitted_strata import fitted_stratum_data
from conditional_remainders import (reuse_source_design_for_evaluation,conditional_design_coefficients,
    conditional_remainder_budgets,full_target_stratified_interval)
from training_average import binomial_ratio_mgf
from joint_contrast import infer_joint_contrast
from repairs import infer_repair


def make_plan(confirmation=False):
    patterns = ['weak_boundary','strong'] if confirmation else ['all_valid','weak_boundary','strong','weak_disjoint','arm_cancel']
    counts = [128,512,2048] if confirmation else [32,128,512]
    return [dict(cell_id=i+1,source_count=count,pattern=pattern,shared_scale=.3,
        n_per_site=1000,source_training_size=500,source_evaluation_size=500,
        target_training_size=500,target_evaluation_size=500,source_design_evaluation_overlap=500,
        valid_fraction=.75,seed_offset=14130000 if confirmation else 13130000)
        for i,(count,pattern) in enumerate((count,pattern) for count in counts for pattern in patterns)]


def run_cell(task):
    setting,repeats,output,envelope = task
    started = time.time()
    data = reuse_source_design_for_evaluation(fitted_stratum_data(setting,repeats))
    valid = int(math.ceil(.75*setting['source_count']))
    coefficients = conditional_design_coefficients(data,valid,.001,.002)
    results = {}
    for noise in ['separate','joint']:
        for bias in ['separate','joint']:
            name = noise+'_noise_'+bias+'_bias'
            results[name] = infer_joint_contrast(data,valid,noise_mode=noise,nuisance_mode=bias,coefficients=coefficients)
    results['joint_mgf_p004'] = infer_joint_contrast(data,valid,source_calibration='training_mgf',envelope=envelope)
    old_budget = conditional_remainder_budgets(data,valid)
    results['previous_arm_projection'] = infer_repair(data,target_calibration='empirical',target_alpha=.02,
        bounded=True,bias_failure=.005,variance_failure=.005,source_mode='dispersion',
        dispersion_calibration='bounded_exponential',source_bias_budget=old_budget['source_mean_budget'],
        target_bias_budget=old_budget['target_budget'])
    reference = full_target_stratified_interval(data,.05)
    reference_width = np.diff(reference['interval'],axis=1)[:,0]
    conditional_mean = data['target_truth'][:,0]+np.mean(data['source_expected'][:,:,0]-data['source_expected'][:,:,1],axis=1)
    actual_drift = np.zeros(repeats)
    for arm,sign in [(0,1),(1,-1)]:
        fixed_valid = np.flatnonzero(data['structural_shift'][:,arm]==0)[:valid]
        actual_drift += sign*data['nuisance_drift'][:,fixed_valid,arm].mean(axis=1)
    actual_dispersion = np.var(data['source_expected'],axis=1)
    rows=[]
    for name,result in results.items():
        interval = result['interval'];width=np.diff(interval,axis=1)[:,0]
        covers = (interval[:,0]<=.1)&(.1<=interval[:,1])
        joint_method = 'source_interval' in result
        source_interval = result.get('source_interval',np.full_like(interval,np.nan))
        source_width = np.diff(source_interval,axis=1)[:,0]
        row=dict(setting,method=name,repeats=repeats,coverage=float(covers.mean()),
            mean_length=float(width.mean()),length_mcse=float(width.std(ddof=1)/np.sqrt(repeats)),
            source_length=float(source_width.mean()),
            source_coverage=float(np.mean((source_interval[:,0]<=.1)&(.1<=source_interval[:,1]))) if joint_method else None,
            rmse=float(np.sqrt(np.mean((result['point_midpoint']-.1)**2))),bias=float(np.mean(result['point_midpoint']-.1)),
            target_length=float(reference_width.mean()),length_ratio_target=float(width.mean()/reference_width.mean()),
            target_rmse=float(np.sqrt(np.mean((reference['point']-.1)**2))),
            target_coverage=float(np.mean((reference['interval'][:,0]<=.1)&(.1<=reference['interval'][:,1]))),
            fallback=float(np.mean(result['fallback']!=0)),
            noise_radius=float(result['noise_radius'].mean()) if joint_method else None,
            nuisance_radius=float(result['nuisance_radius'].mean()) if joint_method else None,
            dispersion_radius=float(result['dispersion_radius'].mean()) if joint_method else None,
            noise_failure=float(np.mean(abs(result['source_center']-conditional_mean)>result['noise_radius']+1e-12)) if joint_method else None,
            nuisance_failure=float(np.mean(abs(actual_drift)>result['nuisance_radius']+1e-12)) if joint_method else None,
            dispersion_failure=float(np.mean(np.any(actual_dispersion>result['dispersion_upper']+1e-12,axis=1))) if joint_method else None)
        rows.append(row)
    names=list(results)
    path=Path(output)/'cells'/('cell_%03d'%setting['cell_id'])
    np.savez_compressed(str(path)+'.npz',methods=np.array(names),
        interval=np.stack([r['interval'] for r in results.values()],axis=1),
        source_interval=np.stack([r.get('source_interval',np.full((repeats,2),np.nan)) for r in results.values()],axis=1),
        midpoint=np.column_stack([r['point_midpoint'] for r in results.values()]),
        source_point=np.column_stack([r['point_source'] for r in results.values()]),
        fallback=np.column_stack([r['fallback'] for r in results.values()]),
        noise_radius=np.column_stack([r.get('noise_radius',np.full(repeats,np.nan)) for r in results.values()]),
        nuisance_radius=np.column_stack([r.get('nuisance_radius',np.full(repeats,np.nan)) for r in results.values()]),
        dispersion_radius=np.column_stack([r.get('dispersion_radius',np.full(repeats,np.nan)) for r in results.values()]),
        source_center=data['target'][:,0]+np.mean(data['source'][:,:,0]-data['source'][:,:,1],axis=1),
        conditional_mean=conditional_mean,actual_nuisance_contrast=actual_drift,
        target_reference_interval=reference['interval'],target_reference_point=reference['point'])
    path.with_suffix('.json').write_text(json.dumps(dict(setting=setting,metrics=rows,runtime_seconds=time.time()-started),indent=2)+'\n')
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output')
    parser.add_argument('--repeats',type=int,default=300)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--confirmation',action='store_true')
    parser.add_argument('--smoke',action='store_true')
    args=parser.parse_args();output=Path(args.output)
    if output.exists() or not str(output).startswith('/scratch.global/zhan9381/FACE-HD/'):
        raise ValueError('Use a new FACE-HD scratch directory')
    if args.repeats<2 or not 1<=args.workers<=4:
        raise ValueError('At least two repetitions and one to four workers required')
    output.mkdir(parents=True);(output/'cells').mkdir()
    plan=make_plan(args.confirmation)
    if args.smoke:plan=plan[:1]
    configuration=dict(settings=plan,repeats=args.repeats,
        confidence_budget=dict(noise=.025,dispersion=.01,nuisance=.005,target=.01),
        scope='Finite strata and design-only weights; arm-specific validity maintained; no high-dimensional claim.',
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')})
    (output/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    envelope=binomial_ratio_mgf(500,4,.04)
    np.savez_compressed(str(output/'training_mgf_envelope.npz'),**envelope)
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for rows in pool.imap_unordered(run_cell,[(s,args.repeats,str(output),envelope) for s in plan]):
            records.extend(rows);print('Completed joint cell',rows[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader()
        writer.writerows(sorted(records,key=lambda r:(r['cell_id'],r['method'])))
    (output/'COMPLETE').write_text('Every source/report interval and outcome retained.\n')


if __name__=='__main__':main()

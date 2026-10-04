"""Paired contrast-dispersion extension and range-adaptive target audit.

These are follow-up analyses of the already saved datasets, not additional
Monte Carlo repetitions or a new independent confirmation panel.
"""
import argparse
import csv
import json
import multiprocessing
from pathlib import Path
import numpy as np
from fitted_strata import fitted_stratum_data
from conditional_remainders import reuse_source_design_for_evaluation,full_target_stratified_interval
from joint_contrast import infer_joint_contrast


def run_cell(task):
    root,profile,setting,repeats=task
    root=Path(root)
    data=reuse_source_design_for_evaluation(fitted_stratum_data(setting,repeats))
    count=setting['source_count'];valid=int(np.ceil(.75*count));overlap=2*valid-count
    previous=np.load(str(root/profile/'cells'/('cell_%03d.npz'%setting['cell_id'])))
    control=infer_joint_contrast(data,valid)
    index=list(previous['methods']).index('joint_noise_joint_bias')
    np.testing.assert_allclose(control['interval'],previous['interval'][:,index],atol=1e-12,rtol=1e-12)
    envelope=dict(np.load(str(root/profile/'training_mgf_envelope.npz')))
    target=full_target_stratified_interval(data,.05,range_adaptive=True)
    results={}
    results['contrast_dispersion_design']=infer_joint_contrast(data,valid,dispersion_mode='contrast')
    results['contrast_dispersion_mgf_p004']=infer_joint_contrast(data,valid,dispersion_mode='contrast',
        source_calibration='training_mgf',envelope=envelope)
    common_valid=np.flatnonzero(np.all(data['structural_shift']==0,axis=1))[:overlap]
    actual_nuisance=np.mean(data['nuisance_drift'][:,common_valid,0]-data['nuisance_drift'][:,common_valid,1],axis=1)
    actual_dispersion=np.var(data['source_expected'][:,:,0]-data['source_expected'][:,:,1],axis=1)
    rows=[]
    for name,result in results.items():
        width=np.diff(result['interval'],axis=1)[:,0]
        source_width=np.diff(result['source_interval'],axis=1)[:,0]
        rows.append(dict(profile=profile,**setting,method=name,repeats=repeats,
            coverage=float(np.mean((result['interval'][:,0]<=.1)&(.1<=result['interval'][:,1]))),
            mean_length=float(width.mean()),source_length=float(source_width.mean()),
            rmse=float(np.sqrt(np.mean((result['point_midpoint']-.1)**2))),bias=float(np.mean(result['point_midpoint']-.1)),
            target_range_length=float(np.diff(target['interval'],axis=1).mean()),
            target_range_coverage=float(np.mean((target['interval'][:,0]<=.1)&(.1<=target['interval'][:,1]))),
            estimated_effect_range=float(target['contrast_range'].mean()),
            noise_radius=float(result['noise_radius'].mean()),nuisance_radius=float(result['nuisance_radius'].mean()),
            dispersion_radius=float(result['dispersion_radius'].mean()),
            nuisance_failure=float(np.mean(abs(actual_nuisance)>result['nuisance_radius']+1e-12)),
            dispersion_failure=float(np.mean(actual_dispersion>result['dispersion_upper'][:,0]+1e-12))))
    output=root/'followup'/profile
    path=output/('cell_%03d'%setting['cell_id'])
    np.savez_compressed(str(path)+'.npz',methods=np.array(list(results)),
        interval=np.stack([r['interval'] for r in results.values()],axis=1),
        source_interval=np.stack([r['source_interval'] for r in results.values()],axis=1),
        midpoint=np.column_stack([r['point_midpoint'] for r in results.values()]),
        noise_radius=np.column_stack([r['noise_radius'] for r in results.values()]),
        nuisance_radius=np.column_stack([r['nuisance_radius'] for r in results.values()]),
        dispersion_radius=np.column_stack([r['dispersion_radius'] for r in results.values()]),
        target_range_interval=target['interval'],target_range_point=target['point'],
        contrast_range=target['contrast_range'],actual_nuisance_contrast=actual_nuisance)
    path.with_suffix('.json').write_text(json.dumps(dict(setting=setting,metrics=rows),indent=2)+'\n')
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('root');parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args();root=Path(args.root)
    if not str(root).startswith('/scratch.global/zhan9381/FACE-HD/') or not 1<=args.workers<=4:
        raise ValueError('Use FACE-HD scratch and one to four workers')
    output=root/'followup';output.mkdir()
    tasks=[]
    for profile in ['main','confirmation']:
        (output/profile).mkdir()
        config=json.loads((root/profile/'configuration.json').read_text())
        tasks.extend((str(root),profile,s,config['repeats']) for s in config['settings'])
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for rows in pool.imap_unordered(run_cell,tasks):
            records.extend(rows);print('Follow-up',rows[0]['profile'],rows[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader()
        writer.writerows(sorted(records,key=lambda r:(r['profile'],r['cell_id'],r['method'])))
    (output/'COMPLETE').write_text('Paired reuse only. Original joint-direction outputs reproduced within 1e-12 for every setting.\n')


if __name__=='__main__':main()

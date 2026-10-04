"""Exact conditional balance diagnostic on paired data and fresh confirmations."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time
import numpy as np
from fitted_strata import fitted_stratum_data
from conditional_standardization import conditional_standardization_interval
from conditional_remainders import full_target_stratified_interval


def run_cell(task):
    root,profile,setting,repeats,output=task
    started=time.time()
    data=fitted_stratum_data(setting,repeats,retain_source_counts=True)
    if profile not in ['fresh','point_check','weak_stress']:
        previous=np.load(str(Path(root)/profile/'cells'/('cell_%03d.npz'%setting['cell_id'])))
        point=full_target_stratified_interval(data)['point']
        np.testing.assert_allclose(point,previous['target_reference_point'],atol=1e-12,rtol=1e-12)
    valid=int(np.ceil(.75*setting['source_count']))
    results={}
    for cohort in ['heldout','all']:
        for mode in (['adaptive'] if profile in ['point_check','weak_stress'] else ['arms','contrast','adaptive']):
            results[cohort+'_'+mode]=conditional_standardization_interval(data,valid,cohort,mode)
    rows=[]
    def diagnostic_mean(values, available):
        return float(np.mean(values[available])) if np.any(available) else None
    for name,result in results.items():
        interval=result['interval'];target=result['target_interval']
        width=np.diff(interval,axis=1)[:,0];target_width=np.diff(target,axis=1)[:,0]
        rows.append(dict(profile=profile,**setting,method=name,repeats=repeats,
            coverage=float(np.mean((interval[:,0]<=.1)&(.1<=interval[:,1]))),
            naive_coverage=float(np.mean((result['naive_interval'][:,0]<=.1)&(.1<=result['naive_interval'][:,1]))),
            naive_length=float(np.diff(result['naive_interval'],axis=1).mean()),
            mean_length=float(width.mean()),length_mcse=float(width.std(ddof=1)/np.sqrt(repeats)),
            rmse=float(np.sqrt(np.mean((result['point_midpoint']-.1)**2))),
            bias=float(np.mean(result['point_midpoint']-.1)),
            projection_rmse=float(np.sqrt(np.mean((result['point_target_projection']-.1)**2))),
            projection_bias=float(np.mean(result['point_target_projection']-.1)),
            source_center_rmse=float(np.sqrt(np.nanmean((result['source_center']-.1)**2))) if np.any(~result['unavailable']) else None,
            source_center_bias=diagnostic_mean(result['source_center']-.1,~result['unavailable']),
            target_length=float(target_width.mean()),length_ratio_target=float(width.mean()/target_width.mean()),
            target_coverage=float(np.mean((target[:,0]<=.1)&(.1<=target[:,1]))),
            target_rmse=float(np.sqrt(np.mean((result['target_point']-.1)**2))),
            noise_radius=diagnostic_mean(result['noise_radius'],~result['unavailable']),
            dispersion_radius=diagnostic_mean(result['dispersion_radius'],~result['unavailable']),
            composition_radius=float(result['composition_radius'].mean()),
            unavailable=float(result['unavailable'].mean()),fallback=float(result['fallback'].mean())))
    path=Path(output)/profile/('cell_%03d'%setting['cell_id'])
    np.savez_compressed(str(path)+'.npz',methods=np.array(list(results)),
        interval=np.stack([r['interval'] for r in results.values()],axis=1),
        midpoint=np.column_stack([r['point_midpoint'] for r in results.values()]),
        target_projection=np.column_stack([r['point_target_projection'] for r in results.values()]),
        naive_interval=np.stack([r['naive_interval'] for r in results.values()],axis=1),
        source_center=np.column_stack([r['source_center'] for r in results.values()]),
        noise_radius=np.column_stack([r['noise_radius'] for r in results.values()]),
        dispersion_radius=np.column_stack([r['dispersion_radius'] for r in results.values()]),
        composition_radius=np.column_stack([r['composition_radius'] for r in results.values()]),
        unavailable=np.column_stack([r['unavailable'] for r in results.values()]),
        fallback=np.column_stack([r['fallback'] for r in results.values()]),
        target_interval=next(iter(results.values()))['target_interval'],target_point=next(iter(results.values()))['target_point'])
    path.with_suffix('.json').write_text(json.dumps(dict(setting=setting,metrics=rows,runtime_seconds=time.time()-started),indent=2)+'\n')
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('root');parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--fresh',action='store_true');parser.add_argument('--point-check',action='store_true')
    parser.add_argument('--weak-stress',action='store_true')
    parser.add_argument('--repeats',type=int,default=500)
    args=parser.parse_args();root=Path(args.root)
    if not str(root).startswith('/scratch.global/zhan9381/FACE-HD/') or not 1<=args.workers<=4 or args.repeats<2:
        raise ValueError('Use FACE-HD scratch and one to four workers')
    output=root/('conditional_weak_stress' if args.weak_stress else ('conditional_point_fresh' if args.point_check else ('conditional_balance_fresh' if args.fresh else 'conditional_balance_paired')))
    output.mkdir()
    tasks=[]
    if args.weak_stress:
        if args.point_check or args.fresh:raise ValueError('Choose exactly one fresh study profile')
        plan=[dict(cell_id=i+1,source_count=count,pattern=pattern,shared_scale=0.,n_per_site=1000,
                   seed_offset=18130000,valid_fraction=.75)
              for i,(count,pattern) in enumerate((count,pattern) for count in [128,512,2048]
                  for pattern in ['all_valid','weak_boundary'])]
        (output/'weak_stress').mkdir()
        tasks=[(str(root),'weak_stress',s,args.repeats,str(output)) for s in plan]
    elif args.point_check:
        if args.fresh:raise ValueError('Choose the score confirmation or the point-rule confirmation')
        plan=[dict(cell_id=i+1,source_count=512,pattern=pattern,shared_scale=.3,n_per_site=1000,
                   seed_offset=17130000,valid_fraction=.75)
              for i,pattern in enumerate(['weak_boundary','strong','arm_cancel'])]
        (output/'point_check').mkdir()
        tasks=[(str(root),'point_check',s,args.repeats,str(output)) for s in plan]
    elif args.fresh:
        plan=[dict(cell_id=i+1,source_count=count,pattern=pattern,shared_scale=.3,n_per_site=1000,
                   seed_offset=16130000,valid_fraction=.75)
              for i,(count,pattern) in enumerate([(128,'weak_boundary'),(512,'weak_boundary'),
                  (2048,'weak_boundary'),(512,'strong'),(512,'arm_cancel')])]
        (output/'fresh').mkdir()
        tasks=[(str(root),'fresh',s,args.repeats,str(output)) for s in plan]
    else:
        plan=[]
        for profile in ['main','confirmation']:
            (output/profile).mkdir()
            config=json.loads((root/profile/'configuration.json').read_text())
            plan.extend(dict(s,profile=profile,repeats=config['repeats']) for s in config['settings'])
            tasks.extend((str(root),profile,s,config['repeats'],str(output)) for s in config['settings'])
    (output/'configuration.json').write_text(json.dumps(dict(settings=plan,fresh=args.fresh or args.point_check or args.weak_stress,point_check=args.point_check,weak_stress=args.weak_stress,
        repeats=args.repeats if args.fresh or args.point_check or args.weak_stress else None,
        confidence_budget=dict(range=.005,composition=.01,source_noise=.01,target_noise=.01,dispersion=.015),
        scope='Exploratory saturated exact-design calibration; valid sites share target conditional means; no high-dimensional model-robustness claim.',
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}),indent=2)+'\n')
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for rows in pool.imap_unordered(run_cell,tasks):
            records.extend(rows);print('Conditional balance',rows[0]['profile'],rows[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted(set(key for row in records for key in row)))
        writer.writeheader();writer.writerows(sorted(records,key=lambda r:(r['profile'],r['cell_id'],r['method'])))
    (output/'COMPLETE').write_text('Every sample and design-fallback case retained.\n')


if __name__=='__main__':main()

"""Paired reanalysis of saved v21 settings under safe budget and design rules."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from scipy.stats import norm

from continuous_balance import UnsupportedDesign
from design_selection import select_source_designs, restrict_weight_inflation, evaluate_retained_sources
from design_budget import apply_design_precision_gate, population_selected_interval
from population_protection import target_coefficient_certificate, composition_radius, expand_conditional_interval
from run_continuous_balance import draw_site
from run_design_followup import draw_design_site, plan
from run_target_protection import fit_packet


VARIANTS = [('legacy', 'legacy', None, False),
            ('conditional_reallocation', 'conditional_reallocation', None, False),
            ('shared_composition', 'shared_composition', None, False),
            ('guard2', 'shared_composition', 2., False),
            ('guard4', 'shared_composition', 4., False),
            ('guard8', 'shared_composition', 8., False),
            ('precision_gate', 'shared_composition', None, True),
            ('guard4_gate', 'shared_composition', 4., True)]
METHODS = [v[0] for v in VARIANTS]+['target_protected', 'target_normal']


def run_repeat(setting, repeat):
    confirmation = setting['panel'] == 'confirmation'
    offset = setting['seed_offset'] if confirmation else 31190000
    rng = np.random.RandomState(offset+10000*setting['cell_id']+repeat)
    draw = draw_site if confirmation else draw_design_site
    target_site = draw(rng, setting, target=True)
    features = target_site[0]
    mean = np.r_[1., features.mean(axis=0)]
    truth = float(np.mean(target_site[3][:, 0]-target_site[3][:, 1]))
    try:
        target, _, _ = fit_packet(target_site, mean)
        certificate = target_coefficient_certificate(features, target_site[1], target_site[2])
    except UnsupportedDesign:
        return dict(methods={name: dict(interval=[-1.,1.], conditional=[-1.,1.], point=0.,
            design_fallback=True, rejected_count=setting['source_count'], gate_rejected=False,
            source_noise=0., dispersion=0., composition=1.) for name in METHODS},
            conditional_truth=truth, target_failed=True)
    outcomes = []
    def source_designs():
        for j in range(setting['source_count']):
            site = draw(rng, setting, index=j)
            outcomes.append(site[2])
            yield np.column_stack((np.ones(len(site[0])), site[0])), site[1]
    base = select_source_designs(source_designs(), mean, 3*setting['source_count']//4,
        design_policy='remove')
    source = evaluate_retained_sources(base, outcomes)
    selections = {limit: restrict_weight_inflation(base, limit) for limit in [None, 2., 4., 8.]}
    results = {}
    for name, budget_policy, limit, use_gate in VARIANTS:
        selected = selections[limit]
        if use_gate:
            selected = apply_design_precision_gate(selected, target['weight_square_sum'],
                model_mode=setting['model_mode'], budget_policy=budget_policy)
        positions = np.searchsorted(base['retained_indices'], selected['retained_indices'])
        subset = None if source is None else {key: value[positions] for key, value in source.items()}
        result = population_selected_interval(selected, subset, target, certificate, len(features),
            model_mode=setting['model_mode'], budget_policy=budget_policy)
        results[name] = dict(interval=result['interval'], conditional=result['conditional_interval'],
            point=result['point'], design_fallback=result['design_fallback'],
            rejected_count=len(selected['rejected']), gate_rejected=selected.get('precision_gate_rejected',False),
            source_noise=result['source_noise'], dispersion=result['dispersion_radius'], composition=result['composition_radius'])
    point = float(target['estimate'][0]-target['estimate'][1])
    noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/.0225))
    conditional = point+noise*np.array([-1.,1.])
    radius = composition_radius(certificate, len(features), .0225)
    results['target_protected'] = dict(interval=expand_conditional_interval(conditional,radius),
        conditional=conditional,point=point,design_fallback=False,rejected_count=0,
        gate_rejected=False,source_noise=0.,dispersion=0.,composition=radius)
    variance = max(0.,target['variance'].sum()+np.var(features.dot(certificate['coefficient']),ddof=1)/len(features))
    results['target_normal'] = dict(interval=np.clip(point+norm.isf(.025)*np.sqrt(variance)*np.array([-1.,1.]),-1,1),
        conditional=point+norm.isf(.025)*np.sqrt(max(0.,target['variance'].sum()))*np.array([-1.,1.]),
        point=point,design_fallback=False,rejected_count=0,gate_rejected=False,source_noise=0.,dispersion=0.,composition=0.)
    return dict(methods=results,conditional_truth=truth,target_failed=False)


def run_cell(task):
    setting,repeats,output,previous = task
    started = time.time()
    records = [run_repeat(setting,r) for r in range(repeats)]
    arrays = {key:np.array([[r['methods'][m][key] for m in METHODS] for r in records])
        for key in ['interval','conditional','point','design_fallback','rejected_count','gate_rejected','source_noise','dispersion','composition']}
    arrays['conditional_truth'] = np.array([r['conditional_truth'] for r in records])
    if setting['model_mode']=='linear':
        oldstem=Path(previous)/setting['panel']/('cell_%03d'%setting['cell_id'])
        old=json.loads(oldstem.with_suffix('.json').read_text());saved=np.load(str(oldstem.with_suffix('.npz')))
        name='source_bernstein' if setting['panel']=='confirmation' else 'remove_linear'
        j=old['methods'].index(name)
        for key in ['interval','conditional','point']:
            np.testing.assert_array_equal(arrays[key][:,0],saved[key][:repeats,j])
    # Every design-fallback shared-budget interval equals the full target comparator.
    for method in ['shared_composition','guard2','guard4','guard8','precision_gate','guard4_gate']:
        j=METHODS.index(method); mask=arrays['design_fallback'][:,j]
        np.testing.assert_array_equal(arrays['interval'][mask,j],arrays['interval'][mask,METHODS.index('target_protected')])
    rows=[]
    for j,name in enumerate(METHODS):
        interval=arrays['interval'][:,j];conditional=arrays['conditional'][:,j]
        rows.append(dict(setting,method=name,repeats=repeats,
            coverage=float(np.mean((interval[:,0]<=.1)&(.1<=interval[:,1]))),
            conditional_coverage=float(np.mean((conditional[:,0]<=arrays['conditional_truth'])&(arrays['conditional_truth']<=conditional[:,1]))),
            length=float(np.diff(interval,axis=1).mean()),rmse=float(np.sqrt(np.mean((arrays['point'][:,j]-.1)**2))),
            fallback_rate=float(arrays['design_fallback'][:,j].mean()),gate_rate=float(arrays['gate_rejected'][:,j].mean()),
            mean_removed=float(arrays['rejected_count'][:,j].mean())))
    stem=Path(output)/(setting['panel']+'_%03d'%setting['cell_id'])
    np.savez_compressed(str(stem)+'.npz',**arrays)
    stem.with_suffix('.json').write_text(json.dumps(dict(setting=setting,methods=METHODS,metrics=rows,
        runtime_seconds=time.time()-started,target_failures=sum(r['target_failed'] for r in records)),indent=2)+'\n')
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--profile',choices=['removal','confirmation','all'],default='all')
    parser.add_argument('--model-mode',choices=['linear','bounded_invalid'],default='linear')
    parser.add_argument('--repeats',type=int,default=200);parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--run-name',default='budget')
    args=parser.parse_args()
    if (Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents
            or not 2<=args.repeats<=200 or not 1<=args.workers<=4
            or Path(args.run_name).name!=args.run_name or args.run_name in ['.','..']):
        raise ValueError('Use a new FACE-HD scratch directory, 2--200 paired repeats and 1--4 workers')
    output=args.root/args.run_name;output.mkdir()
    panels=['confirmation','removal'] if args.profile=='all' else [args.profile]
    settings=[dict(s,panel=panel,model_mode=args.model_mode) for panel in panels for s in plan(panel)]
    (output/'configuration.json').write_text(json.dumps(dict(settings=settings,repeats=args.repeats,
        variants=VARIANTS,previous=str(args.previous),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')}),indent=2)+'\n')
    rows=[]
    with multiprocessing.Pool(args.workers) as pool:
        for result in pool.imap_unordered(run_cell,[(s,args.repeats,str(output),str(args.previous)) for s in settings]):
            rows.extend(result);print('Budget panel',result[0]['panel'],result[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted({key for row in rows for key in row}))
        writer.writeheader();writer.writerows(rows)
    (output/'COMPLETE').write_text('Legacy intervals replay exactly; shared-budget design fallbacks equal protected target.\n')


if __name__=='__main__':main()

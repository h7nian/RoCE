"""Paired target-protection and controlled-misspecification experiments."""
import argparse
import csv
import hashlib
import json
import itertools
import multiprocessing
from pathlib import Path
import time

import numpy as np
from scipy.stats import norm

from continuous_balance import (UnsupportedDesign, linear_balance_design, evaluate_linear_balance,
    conditional_linear_bands, continuous_balance_interval, leaveout_dispersion_upper)
from population_protection import (target_coefficient_certificate, composition_radius,
    expand_conditional_interval, projection_geometry, projection_error_budget)
from run_continuous_balance import draw_site, SUMMARY_FIELDS


METHODS = ['source_universal', 'source_range', 'source_bernstein', 'source_adaptive',
           'target_universal', 'target_range', 'target_bernstein', 'target_adaptive',
           'source_legacy', 'target_legacy', 'target_normal']


def fit_packet(site, target_mean, include_projection=False):
    features, treatment, outcomes, probabilities = site
    matrix = np.column_stack((np.ones(len(features)), features))
    summaries, geometry, diagnostics = [], [], []
    for arm in [1, 0]:
        selected = treatment == arm
        design = linear_balance_design(matrix[selected], target_mean)
        summaries.append(evaluate_linear_balance(design, outcomes[selected]))
        detail = dict(sample_size=int(selected.sum()))
        if include_projection:
            detail.update(projection_geometry(design))
        geometry.append(detail)
        means = probabilities[selected, 1-arm]
        centered = means-.5
        residual = centered-design['basis'].dot(design['basis'].T.dot(centered))
        diagnostics.append(dict(mean=design['weights'].dot(means),
            variance=np.sum(design['weights']**2*means*(1-means)),
            variance_bias=np.sum(design['diagonal']*centered*residual),
            projection_norm=np.linalg.norm(residual)))
    summary = {name: np.array([entry[name] for entry in summaries]) for name in SUMMARY_FIELDS}
    return summary, geometry, diagnostics


def nonlinear_source(rng, setting, index):
    site = draw_site(rng, setting, index=index)
    features, treatment, _, probabilities = site
    if index < setting['source_count']//4:
        probabilities[:, 0] += setting['nonlinearity']*features[:, 0]*features[:, 1]
    if np.any((probabilities < 0) | (probabilities > 1)):
        raise ArithmeticError('Nonlinear DGP probabilities left [0,1]')
    outcomes = rng.binomial(1, probabilities[np.arange(len(features)), 1-treatment])
    return features, treatment, outcomes, probabilities


def validate_nonlinear_support(setting):
    """Check every signal-support corner through the actual DGP implementation."""
    class SupportCorners:
        def uniform(self, lower, upper, shape):
            corners = np.array(list(itertools.product([-1., 1.], repeat=4)))
            features = np.zeros(shape)
            features[:, :4] = corners[np.arange(shape[0]) % len(corners)]
            return features
        def binomial(self, n, probabilities):
            probabilities = np.asarray(probabilities)
            if np.any((probabilities < 0) | (probabilities > 1)):
                raise ValueError('Outcome/treatment probability violates support at a design corner')
            return np.zeros(probabilities.shape, dtype=int)
    nonlinear_source(SupportCorners(), setting, index=0)


def run_repeat(setting, repeat):
    rng = np.random.RandomState(setting['seed_offset']+10000*setting.get('seed_group',setting['cell_id'])+repeat)
    target_site = draw_site(rng, setting, target=True)
    features = target_site[0]
    target_mean = np.r_[1., features.mean(axis=0)]
    conditional_truth = float(np.mean(target_site[3][:, 0]-target_site[3][:, 1]))
    names = METHODS + (['repaired_bernstein'] if 'nonlinearity' in setting else [])
    if setting.get('bounded_invalid',False):names += ['bounded_invalid_bernstein']
    try:
        target, _, _ = fit_packet(target_site, target_mean)
        certificate = target_coefficient_certificate(features, target_site[1], target_site[2])
    except UnsupportedDesign:
        values = {name: dict(interval=np.array([-1., 1.]), conditional=np.array([-1., 1.]),
                            point=0., source_band=np.array([-1., 1.]), source_noise=0.,
                            dispersion=0., composition=1.) for name in names}
        return dict(methods=values, diagnostics=dict(conditional_truth=conditional_truth,
            target_failed=True, source_failed=True, range_upper=2., variance_upper=1.,
            certificate_failed=False))
    values, geometries, diagnostics = [], [], []
    source_failed = False
    for index in range(setting['source_count']):
        site = nonlinear_source(rng, setting, index) if 'nonlinearity' in setting else draw_site(rng, setting, index=index)
        try:
            value, geometry, audit = fit_packet(site, target_mean, include_projection='nonlinearity' in setting or setting.get('bounded_invalid',False))
        except UnsupportedDesign:
            source_failed = True
            continue
        values.append(value); geometries.append(geometry); diagnostics.append(audit)
    count = setting['source_count']
    target_point = float(target['estimate'][0]-target['estimate'][1])
    target_noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(2/.0225))
    target_conditional = target_point+target_noise*np.array([-1., 1.])
    if source_failed:
        conditional = target_point+np.sqrt(target['weight_square_sum'].sum()/2*np.log(200.))*np.array([-1., 1.])
        bands = dict(conditional_interval=conditional, source_band=np.array([-1., 1.]),
                     source_noise=0., dispersion_radius=0., target_projection=target_point)
    else:
        source = {key: np.stack([v[key] for v in values]) for key in SUMMARY_FIELDS}
        bands = conditional_linear_bands(source, target, 3*count//4, dispersion_failure=.015)
    methods = {}
    def add(name, conditional, radius, point, source_band, noise=0., dispersion=0.):
        interval = expand_conditional_interval(conditional, radius)
        methods[name] = dict(interval=interval, conditional=conditional, point=float(np.clip(point, interval[0], interval[1])),
            source_band=source_band, source_noise=noise, dispersion=dispersion, composition=radius)
    for mode in ['universal', 'range', 'bernstein', 'adaptive']:
        add('source_'+mode, bands['conditional_interval'], composition_radius(certificate, len(features), .01, mode),
            bands['target_projection'], bands['source_band'], bands['source_noise'], bands['dispersion_radius'])
        add('target_'+mode, target_conditional, composition_radius(certificate, len(features), .0225, mode),
            target_point, target_conditional)
    legacy = continuous_balance_interval(source, target, 3*count//4, len(features)) if not source_failed else None
    if legacy is not None:
        add('source_legacy', legacy['conditional_interval'], legacy['composition_radius'],
            legacy['target_projection'], legacy['source_band'], legacy['source_noise'], legacy['dispersion_radius'])
    else:
        methods['source_legacy'] = methods['source_universal'].copy()
    old_noise = np.sqrt(target['weight_square_sum'].sum()/2*np.log(80.))
    add('target_legacy', target_point+old_noise*np.array([-1., 1.]),
        2*np.sqrt(np.log(80.)/(2*len(features))), target_point, target_conditional)
    variance = max(0., target['variance'].sum()+np.var(features.dot(certificate['coefficient']), ddof=1)/len(features))
    normal = target_point+norm.isf(.025)*np.sqrt(variance)*np.array([-1., 1.])
    methods['target_normal'] = dict(interval=np.clip(normal, -1, 1),
        conditional=target_point+norm.isf(.025)*np.sqrt(max(0., target['variance'].sum()))*np.array([-1., 1.]),
        point=target_point, source_band=np.array([-1., 1.]), source_noise=0., dispersion=0., composition=0.)
    audit = dict(conditional_truth=conditional_truth, target_failed=False, source_failed=source_failed,
        range_upper=certificate['range_upper'], variance_upper=certificate['sample_variance_upper'],
        certificate_failed=bool(certificate['range_upper'] < 2*setting['heterogeneity'] or
            certificate['sample_variance_upper'] < np.var(target_site[3][:, 0]-target_site[3][:, 1], ddof=1)))
    if not source_failed:
        means = np.array([[a['mean'] for a in row] for row in diagnostics])
        variance_bias = np.array([[a['variance_bias'] for a in row] for row in diagnostics])
        audit['max_valid_drift'] = float(np.max(np.abs(means[count//4:]-target_site[3].mean(axis=0))))
        audit['variance_centering_bias'] = float((count-1)/count**2*variance_bias.sum())
        audit['contrast_dispersion'] = float(np.var(means[:, 0]-means[:, 1]))
        protection_profiles = []
        if 'nonlinearity' in setting:
            protection_profiles.append(('repaired_bernstein', [.2, 0.]))
        if setting.get('bounded_invalid', False):
            protection_profiles.append(('bounded_invalid_bernstein', [.5, .5]))
        for label, envelopes in protection_profiles:
            budgets = {}
            for arm in range(2):
                summary = {key: value[:, arm] for key, value in source.items()}
                for key in ['sample_size', 'projection_commutator_norm']:
                    summary[key] = np.array([g[arm][key] for g in geometries])
                rho = envelopes[arm]*np.sqrt(summary['sample_size'])
                budgets[('mu1','mu0')[arm]] = projection_error_budget(summary, rho, count//4)
            budgets['contrast'] = {key: budgets['mu1'][key]+budgets['mu0'][key] for key in budgets['mu1']}
            repaired = conditional_linear_bands(source, target, 3*count//4,
                dispersion_failure=.015, dispersion_budgets=budgets)
            add(label, repaired['conditional_interval'], composition_radius(certificate, len(features), .01),
                repaired['target_projection'], repaired['source_band'], repaired['source_noise'], repaired['dispersion_radius'])
            audit[label+'_centering_allowance'] = float((count-1)/count**2*budgets['contrast']['variance_bias_sum'])
    else:
        if 'nonlinearity' in setting:
            methods['repaired_bernstein'] = methods['source_bernstein'].copy()
        if setting.get('bounded_invalid', False):
            methods['bounded_invalid_bernstein'] = methods['source_bernstein'].copy()
    return dict(methods=methods, diagnostics=audit)


def plan(profile):
    if profile=='pilot':
        cells=[dict(dimension=p, source_count=16, heterogeneity=.15, pattern='weak') for p in [4,10]]
    elif profile=='target':
        cells=[dict(dimension=p, source_count=k, heterogeneity=h, pattern='weak')
               for p in [4,10] for k in [64,256] for h in [0.,.15]]
        cells += [dict(dimension=p, source_count=2048, heterogeneity=0., pattern='weak') for p in [4,10]]
        cells += [dict(dimension=p, source_count=256, heterogeneity=.15, pattern='all_valid') for p in [4,10]]
    elif profile=='misspecification':
        cells=[dict(dimension=4, source_count=k, heterogeneity=.15, pattern='weak', nonlinearity=eta, seed_group=i+1)
               for i,k in enumerate([64,256]) for eta in [0.,.05,.1,.15]]
    else:
        raise ValueError('Unknown profile')
    return [dict(c, cell_id=i+1, n_per_site=1000, source_shift=.4,
                 seed_offset={'pilot':23190000,'target':24190000,'misspecification':25190000}[profile]) for i,c in enumerate(cells)]


def run_cell(task):
    setting, repeats, output = task
    started=time.time(); records=[run_repeat(setting,r) for r in range(repeats)]
    methods=list(records[0]['methods']); diagnostics=[r['diagnostics'] for r in records]
    arrays={name:np.array([[r['methods'][m][name] for m in methods] for r in records])
            for name in ['interval','conditional','point','source_band','source_noise','dispersion','composition']}
    arrays['conditional_truth']=np.array([d['conditional_truth'] for d in diagnostics])
    prefix=Path(output)/('cell_%03d'%setting['cell_id']); np.savez_compressed(str(prefix)+'.npz',**arrays)
    metrics=[]
    for index,method in enumerate(methods):
        interval=arrays['interval'][:,index]; conditional=arrays['conditional'][:,index]
        coverage=(interval[:,0]<=.1)&(.1<=interval[:,1])
        conditional_coverage=(conditional[:,0]<=arrays['conditional_truth'])&(arrays['conditional_truth']<=conditional[:,1])
        errors=arrays['point'][:,index]-.1
        metrics.append(dict(setting,method=method,repeats=repeats,coverage=float(coverage.mean()),
            conditional_coverage=float(conditional_coverage.mean()),length=float(np.diff(interval,axis=1).mean()),
            conditional_length=float(np.diff(conditional,axis=1).mean()),rmse=float(np.sqrt(np.mean(errors**2))),
            bias=float(errors.mean()),composition_radius=float(arrays['composition'][:,index].mean()),
            source_noise=float(arrays['source_noise'][:,index].mean()),dispersion_radius=float(arrays['dispersion'][:,index].mean())))
    Path(str(prefix)+'.json').write_text(json.dumps(dict(setting=setting,methods=methods,metrics=metrics,
        diagnostics=diagnostics,runtime_seconds=time.time()-started),indent=2)+'\n')
    return metrics


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--profile',choices=['pilot','target','misspecification'],default='target')
    parser.add_argument('--repeats',type=int,default=200);parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--run-name',help='New output-directory name; defaults to the profile')
    parser.add_argument('--bounded-invalid',action='store_true',help='Also protect arbitrary bounded invalid means')
    args=parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats<2 or not 1<=args.workers<=4:
        raise ValueError('Use FACE-HD scratch, >=2 repetitions, and 1--4 workers')
    name=args.run_name or args.profile
    if Path(name).name!=name or name in ['.','..']:raise ValueError('Use one output-directory name')
    output=args.root/name;output.mkdir();settings=plan(args.profile)
    for setting in settings:
        if args.bounded_invalid:setting['bounded_invalid']=True
        if 'nonlinearity' in setting:validate_nonlinear_support(setting)
    configuration=dict(settings=settings,repeats=args.repeats,
        ledger=dict(coefficient=.005,composition=.01,source_outcome=.01,target_outcome=.01,dispersion=.015),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')})
    (output/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    rows=[]
    with multiprocessing.Pool(args.workers) as pool:
        for result in pool.imap_unordered(run_cell,[(s,args.repeats,str(output)) for s in settings]):
            rows.extend(result);print('Completed',args.profile,result[0]['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted({k for r in rows for k in r}));writer.writeheader();writer.writerows(rows)
    (output/'COMPLETE').write_text('Every repetition and all conditional/population intervals retained.\n')


if __name__=='__main__':main()

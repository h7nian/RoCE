"""Audit saved exact-balance intervals against conditional and population truth."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import binom, chi2
from fitted_strata import fitted_stratum_data
from conditional_standardization import target_stratified_components


def audit_naive_variance(root):
    cells=np.array([(first,second) for first in [-1,1] for second in [-1,1]])
    means=np.column_stack((.55+.08*cells[:,1],.45+.08*cells[:,1]))
    site_variances=[]
    site_missing_bounds=[]
    for orientation in [-1,1]:
        mass=np.exp(.3*orientation*cells[:,0]);mass/=mass.sum()
        propensity=.5+.1*orientation*cells[:,0]+.05*cells[:,1]
        joint=mass[:,None]*np.column_stack((propensity,1-propensity))
        successes=np.arange(1,1001)
        inverse=np.sum(binom.pmf(successes[:,None,None],1000,joint[None])/successes[:,None,None],axis=0)
        site_variances.append(np.sum((.25**2+.25*.75/1000)*means*(1-means)*inverse))
        site_missing_bounds.append(float(binom.cdf(1,1000,joint).sum()))
    rows=[]
    for cell,count in [(1,128),(3,512),(5,2048)]:
        packet=np.load(str(root/'conditional_weak_stress/weak_stress'/('cell_%03d.npz'%cell)))
        index=list(packet['methods']).index('all_adaptive')
        interval=packet['naive_interval'][:,index]
        if np.any(abs(interval)>=1):raise AssertionError('A clipped radius needs direct variance records')
        estimated=((interval[:,1]-interval[:,0])/(2*1.959963984540054))**2
        empirical=np.var(packet['source_center'][:,index],ddof=1)
        analytic=np.mean(site_variances)/count
        degrees=len(estimated)-1
        rows.append(dict(source_count=count,analytic_variance=analytic,
            mean_estimated_variance=float(estimated.mean()),empirical_variance=float(empirical),
            mean_estimated_over_analytic=float(estimated.mean()/analytic),
            empirical_over_analytic=float(empirical/analytic),
            unavailable_design_probability_upper=count*np.mean(site_missing_bounds),
            normal_reference_two_sided_variance_tail=float(2*min(chi2.cdf(degrees*empirical/analytic,degrees),
                chi2.sf(degrees*empirical/analytic,degrees)))))
    (root/'report/naive_analytic_variance_audit.json').write_text(json.dumps(rows,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('root');root=Path(parser.parse_args().root)
    cache={};rows=[]
    folders=['conditional_balance_paired','conditional_balance_fresh','conditional_point_fresh','conditional_weak_stress']
    for folder in folders:
        configuration=json.loads((root/folder/'configuration.json').read_text())
        for setting in configuration['settings']:
            if folder=='conditional_balance_paired':profile=setting['profile'];repeats=setting['repeats']
            else:
                profile={'conditional_balance_fresh':'fresh','conditional_point_fresh':'point_check','conditional_weak_stress':'weak_stress'}[folder]
                repeats=configuration['repeats']
            key=(setting['seed_offset'],setting['shared_scale'],repeats)
            if key not in cache:
                target_data=fitted_stratum_data(dict(setting,source_count=1,pattern='all_valid'),repeats)
                counts=target_data['target_training_counts']+target_data['target_evaluation_counts']
                masses=counts.sum(axis=(2,3))/counts.sum(axis=(1,2,3))[:,None]
                truth=.1+setting['shared_scale']*(masses@np.array([-1.,-1.,1.,1.]))
                point=target_stratified_components(counts,.005,.0225,.0225)['point']
                cache[key]=(truth,point)
            truth,point=cache[key]
            packet=np.load(str(root/folder/profile/('cell_%03d.npz'%setting['cell_id'])))
            np.testing.assert_allclose(packet['target_point'],point,atol=1e-12)
            count=setting['source_count'];pattern=setting['pattern']
            shift=2*count**(-.25)*np.sqrt(.5/250)
            structural_bias={'all_valid':0.,'weak_boundary':shift/4,'strong':.05,
                             'weak_disjoint':shift/2,'arm_cancel':0.}[pattern]
            for index,method in enumerate(packet['methods']):
                interval=packet['interval'][:,index]
                if np.any(abs(interval)>=1):
                    raise AssertionError('Regenerate clipped cases before reconstructing conditional endpoints')
                population_radius=packet['composition_radius'][:,index]
                cohort_interval=interval-population_radius[:,None]*np.array([-1.,1.])
                projection=np.clip(point,cohort_interval[:,0],cohort_interval[:,1])
                if 'target_projection' in packet:
                    np.testing.assert_allclose(projection,packet['target_projection'][:,index],atol=1e-12)
                source_radius=packet['noise_radius'][:,index]+packet['dispersion_radius'][:,index]
                source_center=packet['source_center'][:,index]
                available=~packet['unavailable'][:,index]
                source_covers=available&(abs(source_center-truth)<=source_radius)
                if np.any(abs(projection[source_covers]-truth[source_covers]) > abs(point[source_covers]-truth[source_covers])+1e-12):
                    raise AssertionError('Target projection violates its conditional contraction property')
                rows.append(dict(folder=folder,profile=profile,cell_id=setting['cell_id'],method=method,
                    source_conditional_coverage=float(np.mean(source_covers)) if available.all() else None,
                    bias_certificate_failure=float(np.mean(packet['dispersion_radius'][:,index]+1e-12<abs(structural_bias))),
                    composition_failure=float(np.mean(abs(truth-.1)>population_radius+1e-12)),
                    projection_rmse=float(np.sqrt(np.mean((projection-.1)**2))),
                    population_coverage=float(np.mean((interval[:,0]<=.1)&(.1<=interval[:,1]))),
                    projection_geometry_failures=0))
    output=root/'report';output.mkdir(exist_ok=True)
    with (output/'conditional_truth_audit.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(method_setting_rows=len(rows),minimum_source_conditional_coverage=min(r['source_conditional_coverage'] for r in rows),
        maximum_bias_certificate_failure=max(r['bias_certificate_failure'] for r in rows),
        maximum_composition_failure=max(r['composition_failure'] for r in rows),projection_geometry_failures=0)
    (output/'conditional_truth_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    audit_naive_variance(root)
    print(json.dumps(summary))


if __name__=='__main__':main()

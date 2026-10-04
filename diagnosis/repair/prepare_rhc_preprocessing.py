#!/usr/bin/env python3
"""Freeze a paired RHC outer-preprocessing comparison on the published folds."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import EXCLUDED_NODES


def prepare(output, checks):
    output, checks = pilot.scratch_path(output), pilot.scratch_path(checks)
    if output.parent != pilot.SCRATCH/'real_data/rhc':
        raise ValueError('Use the RHC scratch hierarchy')
    marker = checks/'RHC_PREPROCESSING_CHECKS_PASSED'
    if not marker.is_file():
        raise ValueError('Outer-preprocessing regression checks must pass first')
    repository = Path(__file__).resolve().parents[2]
    for name, digest in json.loads((checks/'source_manifest.json').read_text()).items():
        if pilot.digest(repository/name) != digest:
            raise ValueError('Code changed after regression checks: '+name)
    library = checks/'Rlib'
    prior_library = pilot.SCRATCH/'implementation/r12/density_cv_retry_v2/Rlib'
    if pilot.digest(library/'RoCE/libs/RoCE.so') != pilot.digest(prior_library/'RoCE/libs/RoCE.so'):
        raise ValueError('This paired comparison must retain the numerical kernel')
    output.mkdir(exist_ok=False)
    for directory in ['workflow','tasks','logs','submissions']:
        (output/directory).mkdir()
    source = Path(__file__).resolve().parent
    for name in ['submit_repeat_pilot.py','run_repeat_pilot.sh','run_rhc_current.R',
                 'layer_comparison_summary.R','rhc_stage_signatures.R','prepare_rhc_preprocessing.py']:
        shutil.copy2(str(source/name),str(output/'workflow'/name))
    shutil.copy2(str(repository/'R/real_data_rhc.R'),str(output/'workflow/real_data_rhc.R'))
    reference = output.parent/'report_selection_private_min_radius3_v1/methods.csv'
    shutil.copy2(str(reference),str(output/'historical_methods.csv'))
    tasks=[]
    for method in ['RoCE']+list(pilot.BASELINE_METHODS):
        for preprocessing in ['cohort','outer_fold']:
            tasks.append(dict(task_id=len(tasks)+1,
                config='RHC_preprocessing_{}_{}'.format(preprocessing,method),K=3,sim_id=42,
                role='method' if method=='RoCE' else 'baseline',
                baseline='' if method=='RoCE' else method,radius=3,target_radius=5,
                nuisance_rule='min',source_program='calibrated',target_program='hou_calibrated',
                fold_seed=0,preprocessing=preprocessing))
    with (output/'manifest.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(tasks[0]));writer.writeheader();writer.writerows(tasks)
    configuration=dict(library=str(library),worker_script='run_rhc_current.R',
        cpus=4,memory='8G',time_limit='01:00:00',checkpoint=True,max_restarts=20,
        account='hou00123',partitions=pilot.DEFAULT_PARTITIONS,exclude_nodes=EXCLUDED_NODES,
        crossfit_layers=2,baseline_methods=list(pilot.BASELINE_METHODS),covariate_profile='historical',
        target_site='Private',source_radius=3,target_radius=5,n_folds=10,nlambda=100,
        aggregation_cutoff=2,n_bootstrap=5000,collation='C.UTF-8',
        data_sha256=pilot.digest(library/'RoCE/extdata/rhc.csv'),
        scope='Paired outer-preprocessing change; fixed cohort, target, fitting rules and original method-specific folds',
        inner_preprocessing_scope='Initial/calibration CV reuse the outer-training transformation; not fully nested preprocessing',
        baseline_scope='SS/IVW exclude their site validation rows; full-sample DR baselines retain their original protocol',
        communication_scope='Federated emulation with pooled covariate preprocessing within allowed rows')
    pilot.write_json(output/'configuration.json',configuration)
    (output/'README.md').write_text('''# Paired RHC outer-preprocessing comparison

Keep all5039 patients,61 working features,Private target,the original partition,
min CV,100 lambdas,Two-layer/one-round,joint TATE,cutoff2,source radius3,target5.
For each method run both whole-cohort and training-only outer preprocessing.
Task1/2 includes RoCE and calibrated Target-only; tasks3--10 cover four baselines.
Each task is a separate Slurm job with checkpoints/requeue. Old runs are untouched.

The new RoCE mode learns a shared fill/center/scale outside all sites' outer k1.
All initial fits,calibration,eta learning andevaluation for that k1 use the
same feature coordinates. Initial k2/CV preprocessing is not further refitted;
this is an outer-boundary sensitivity,not proof of every nested boundary.
The fixed category dictionary is independent of held-out observed categories.

SS/IVW retain their original site-specific cross-fitting RNG streams. Their
preprocessing excludes the current site's validation rows while using other
sites' available covariates. Federated-DR and Pooled-DR fit nuisances on the full
sample,as before,anddo not acquire a new outer cross-fitting procedure here.
Their paired equality is checked. All standard errors retain their original
conditional inference formulas; this experiment does not identify causal bias.

First verify that whole-cohort controls reproduce historical_methods.csv.
Compare all paired point estimates,SEs,CIs,arm means/covariance andsource weights.
No data-specific result is selected or overwritten based on this comparison.
''')
    paths=[output/'manifest.csv',output/'configuration.json',output/'README.md',
           output/'historical_methods.csv',marker,checks/'source_manifest.json']
    paths+=list((output/'workflow').iterdir())
    paths += [p for p in (library/'RoCE').rglob('*') if p.is_file()]
    pilot.write_json(output/'provenance.json',{str(p):pilot.digest(p) for p in paths})
    print('Prepared10 paired tasks; numerical kernel unchanged.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--checks',type=Path,required=True)
    args=parser.parse_args()
    prepare(args.output,args.checks)

"""Random-design source eligibility and matrix-Chernoff diagnostics."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
from scipy.special import expit
from continuous_balance import UnsupportedDesign, linear_balance_design


def chernoff_design_bound(source_count, dimension, source_size, gram_lower_bound,
                          feature_norm_square=None):
    columns=dimension+1
    bound=columns if feature_norm_square is None else feature_norm_square
    if min(source_count,columns,source_size,gram_lower_bound,bound)<=0:
        raise ValueError('Positive design constants are required')
    failure=min(1.,2*source_count*columns*np.exp(-source_size*gram_lower_bound/(8*bound)))
    leverage_bound=2*bound/(source_size*gram_lower_bound)
    return dict(any_source_gram_failure=failure,leverage_upper=leverage_bound,
                condition_upper=np.sqrt(2*bound/gram_lower_bound))


def run_cell(task):
    setting,repeats,output=task;started=time.time();records=[]
    count=setting['source_count'];dimension=setting['dimension'];size=setting['source_size']
    for repeat in range(repeats):
        rng=np.random.RandomState(29190000+setting['cell_id']*10000+repeat)
        target_features=rng.uniform(-1,1,(1000,dimension));target_mean=np.r_[1.,target_features.mean(axis=0)]
        failed=0;norms=[];leverage=[];min_eigenvalues=[]
        for site in range(count):
            features=rng.uniform(-1,1,(size,dimension))
            propensity=expit(setting['overlap_strength']*(features[:,0]+features[:,1]))
            treatment=rng.binomial(1,propensity)
            matrix=np.column_stack((np.ones(size),features))
            site_failed=False
            for arm in [1,0]:
                rows=treatment==arm;selected=matrix[rows]
                gram=selected.T.dot(selected)/size
                min_eigenvalues.append(float(np.linalg.eigvalsh(gram).min()))
                try:
                    design=linear_balance_design(selected,target_mean)
                except UnsupportedDesign:
                    site_failed=True;continue
                norms.append(design['weight_square_sum']);leverage.append(design['leverage'].max())
            failed+=int(site_failed)
        valid=3*count//4
        # This is an eligibility calculation, not an implemented alternative CI.
        adjusted_valid=max(0,valid-failed)
        records.append(dict(failed_sites=failed,any_failure=failed>0,
            valid_after_design_removal=adjusted_valid,retained_sites=count-failed,
            eligible_after_count_adjustment=adjusted_valid>=1,
            max_weight_square_sum=max(norms) if norms else None,
            max_leverage=max(leverage) if leverage else None,min_gram_eigenvalue=min(min_eigenvalues)))
    # Uniform X and propensity >= expit(-2s) give G >= pi_min*diag(1,1/3,...).
    gram_lower=expit(-2*setting['overlap_strength'])/3
    sufficient=chernoff_design_bound(count,dimension,size,gram_lower)
    fractions=np.array([r['failed_sites']/count for r in records])
    any_failure=np.array([r['any_failure'] for r in records])
    row=dict(setting,repeats=repeats,site_failure_rate=float(fractions.mean()),
        site_failure_mcse=float(fractions.std(ddof=1)/np.sqrt(repeats)),
        any_source_failure=float(any_failure.mean()),
        any_failure_mcse=float(np.sqrt(any_failure.mean()*(1-any_failure.mean())/repeats)),
        eligible_after_count_adjustment=float(np.mean([r['eligible_after_count_adjustment'] for r in records])),
        max_observed_leverage=max([r['max_leverage'] for r in records if r['max_leverage'] is not None],default=None),
        gram_lower_assumption=gram_lower,**sufficient)
    path=Path(output)/('cell_%03d.json'%setting['cell_id'])
    path.write_text(json.dumps(dict(setting=setting,metrics=row,records=records,runtime_seconds=time.time()-started),indent=2)+'\n')
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    parser.add_argument('--repeats',type=int,default=200);parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    if Path('/scratch.global/zhan9381/FACE-HD') not in args.root.resolve().parents or args.repeats<2 or not 1<=args.workers<=4:
        raise ValueError('Use FACE-HD scratch, >=2 repeats and 1--4 workers')
    settings=[dict(cell_id=i+1,source_size=n,dimension=p,source_count=k,overlap_strength=s)
              for i,(n,p,k,s) in enumerate((n,p,k,s) for n in [120,200,1000]
                  for p in [10,45] for k in [4,16,64] for s in [0.,2.5])]
    settings += [dict(cell_id=37+i,source_size=n,dimension=4,source_count=k,overlap_strength=0.)
                 for i,(n,k) in enumerate([(1000,4),(4000,16),(16000,64)])]
    output=args.root/'design_stability';output.mkdir()
    (output/'configuration.json').write_text(json.dumps(dict(settings=settings,repeats=args.repeats,
        target_size=1000,seed_offset=29190000,
        scope='Design-only eligibility; count-adjusted removal is reported as a candidate, not an evaluated CI',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    rows=[]
    with multiprocessing.Pool(args.workers) as pool:
        for row in pool.imap_unordered(run_cell,[(s,args.repeats,str(output)) for s in settings]):
            rows.append(row);print('Design stability',row['cell_id'],flush=True)
    with (output/'metrics.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/'COMPLETE').write_text('All design failures retained; no outcome-based filtering.\n')


if __name__=='__main__':main()

"""Reproduce paired target samples and add an unsplit saturated reference."""
import csv
import json
import argparse
from pathlib import Path
import numpy as np
from fitted_strata import fitted_stratum_data
from conditional_remainders import full_target_stratified_interval


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('root')
    root=Path(parser.parse_args().root)
    out=root/'target_baseline_audit';out.mkdir(exist_ok=True)
    cache={};rows=[]
    for profile in ['paired','reused','confirmation','large_k']:
     config_file=root/profile/'configuration.json'
     if not config_file.exists():continue
     config=json.loads(config_file.read_text())
     for setting in config['settings']:
      packet_file=root/profile/'cells'/('cell_%03d.npz'%setting['cell_id'])
      if not packet_file.exists():continue
      key=(setting['seed_offset'],setting['shared_scale'],config['repeats'])
      if key not in cache:
       diagnostic=dict(setting,source_count=1,pattern='all_valid')
       data=fitted_stratum_data(diagnostic,config['repeats'])
       cache[key]=(data,full_target_stratified_interval(data))
      data,result=cache[key]
      packet=np.load(str(packet_file))
      np.testing.assert_allclose(data['reference_point'],packet['target_reference_point'],atol=1e-12)
      interval=result['interval'];width=np.diff(interval,axis=1)[:,0]
      path=out/(profile+'_cell_%03d.npz'%setting['cell_id'])
      np.savez_compressed(str(path),interval=interval,point=result['point'])
      for index,name in enumerate(packet['methods']):
       method_width=np.diff(packet['interval'][:,index],axis=1)[:,0]
       rows.append(dict(profile=profile,cell_id=setting['cell_id'],source_count=setting['source_count'],pattern=setting['pattern'],
        shared_scale=setting['shared_scale'],method=name,target_length=float(width.mean()),
        target_coverage=float(np.mean((interval[:,0]<=.1)&(.1<=interval[:,1]))),
        target_rmse=float(np.sqrt(np.mean((result['point']-.1)**2))),length_ratio=float(method_width.mean()/width.mean())))
    with (out/'metrics.csv').open('w') as f:
     w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('Audited',len(rows),'method-setting rows')


if __name__ == "__main__":
    main()

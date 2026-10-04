#!/usr/bin/env python3
"""Report every paired preprocessing result after all ten tasks complete."""
import argparse
import csv
import json
from pathlib import Path

METHOD_LABELS = {'Calibrated target-only':'Target-only','sample_size':'SS',
                'inverse_variance':'IVW','federated_dr':'Federated-DR',
                'pooled_dr':'Pooled-DR','RoCE':'RoCE'}
FIELDS = ['estimate','se','ci_lower','ci_upper']


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def review(root, output):
    root, output = root.resolve(), output.resolve()
    if root not in output.parents:
        raise ValueError('Keep the review in its RHC study directory')
    historical = {row['method']:row for row in read(root/'historical_methods.csv')}
    results = {}
    control_errors = []
    for task in read(root/'manifest.csv'):
        directory = root/'tasks'/task['task_id']
        if not (directory/'COMPLETE').is_file():
            raise ValueError('Task '+task['task_id']+' is not complete')
        rows = read(directory/'methods.csv')
        wanted = ['Calibrated target-only','RoCE'] if task['role']=='method' else [task['baseline']]
        for method in wanted:
            record = next(row for row in rows if row['method']==method)
            values = {name:float(record[name]) for name in FIELDS}
            if values['se']<=0:
                raise ValueError('Invalid reported standard error')
            if task['preprocessing']=='cohort':
                for name in FIELDS:
                    difference = abs(values[name]-float(historical[method][name]))
                    control_errors.append(difference)
                    if difference>1e-10:
                        raise ValueError('Whole-cohort control differs from historical '+method+': '+name)
            results[method,task['preprocessing']]=values
    paired=[]
    for method,label in METHOD_LABELS.items():
        old,new=results[method,'cohort'],results[method,'outer_fold']
        paired.append(dict(method=label,
          cohort_estimate_pp=100*old['estimate'],outer_estimate_pp=100*new['estimate'],
          change_pp=100*(new['estimate']-old['estimate']),cohort_se_pp=100*old['se'],
          outer_se_pp=100*new['se'],se_change_percent=100*(new['se']/old['se']-1),
          outer_ci_lower_pp=100*new['ci_lower'],outer_ci_upper_pp=100*new['ci_upper']))
    for method in ['federated_dr','pooled_dr']:
        if max(abs(results[method,'cohort'][field]-results[method,'outer_fold'][field]) for field in FIELDS)>1e-10:
            raise ValueError('Full-sample DR baseline changed unexpectedly: '+method)
    output.mkdir(parents=True,exist_ok=False)
    with (output/'paired_methods.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(paired[0]));writer.writeheader();writer.writerows(paired)
    gain=100*(1-results['RoCE','outer_fold']['se']/results['Calibrated target-only','outer_fold']['se'])
    checks=dict(tasks=10,displayed_methods=6,max_control_difference=max(control_errors),
        unchanged_full_sample_dr=True,outer_roce_se_reduction_percent=gain,
        scope='Outer-training preprocessing; inner calibration/CV transforms are not separately refitted')
    (output/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    lines=['# RHC 配对预处理结果','',
      '保留5039例、Private target、原始fold、min、100个lambda、Two-layer、joint TATE、cutoff2、source M=3及target M=5。',
      '旧结果全部保留。本次只比较预处理边界，不估计真实临床bias或coverage。','',
      '| 方法 | 原TATE | 新TATE | 原SE | 新SE | 新95% CI |',
      '|---|---:|---:|---:|---:|---|']
    for r in paired:
        lines.append('| {method} | {cohort_estimate_pp:.3f} | {outer_estimate_pp:.3f} | {cohort_se_pp:.3f} | {outer_se_pp:.3f} | [{outer_ci_lower_pp:.3f}, {outer_ci_upper_pp:.3f}] |'.format(**r))
    lines += ['', '数值单位均为百分点。新RoCE相对新calibrated Target-only的报告SE降低{:.2f}%。'.format(gain),
      '所有全队列对照均复现历史值，最大绝对差异为{:.3g}。'.format(max(control_errors)), '',
      '边界：RoCE和calibrated Target-only使用排除所有site外层k1的共同转换。SS/IVW保持其原本site-specific folds，排除对应site验证行；其他site协变量仍可用于预处理。Federated-DR/Pooled-DR保持原全样本拟合，配对结果应不变。', '',
      '本次并未对每个初始k2或内部CV训练集再次拟合预处理；不能将这一外层边界验证解释为所有嵌套预处理或整个渐近证明已验证。']
    (output/'README_zh.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[4:13]))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();review(args.root,args.output)

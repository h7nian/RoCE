#!/usr/bin/env python3
"""Extend the completed RHC assessment with its clean OR-branch validation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import math
from uuid import uuid4


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def append(base, output):
    base, output = base.resolve(), output.resolve()
    if not str(output).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Use the RHC scratch hierarchy')
    previous = base/'final_assessment_private_min_radius3_v1'
    study = base/'model_validation_or_branch_v1'
    review = study/'reviews/complete_v1'
    audit = study/'reviews/branch_audit_v1'
    figures = study/'reviews/figures_v1'
    for directory in (previous, review, audit, figures):
        if not (directory/'CHECKS_PASSED').is_file():
            raise ValueError('Incomplete review: '+str(directory))
    if not (previous/'VISUAL_REVIEW_PASSED').is_file() or not (figures/'VISUAL_REVIEW_PASSED').is_file():
        raise ValueError('Both the original and supplemental figures require visual review')
    if output.exists():
        raise ValueError('Preserve the existing assessment; use a new output directory')
    summary = read_rows(review/'summary.csv')
    records = read_rows(review/'replicates.csv')
    if len(summary) != 5 or len(records) != 1000:
        raise ValueError('The full200-repeat OR branch is required')
    for row in summary:
        values = [r for r in records if r['label']==row['method']]
        if sorted(int(r['repeat_id']) for r in values) != list(range(1,201)):
            raise ValueError('Missing or repeated seed')
        errors = [float(r['estimate'])-float(r['truth']) for r in values]
        expected = dict(bias=statistics.mean(errors), rmse=math.sqrt(statistics.mean(e*e for e in errors)),
            empirical_sd=statistics.stdev(float(r['estimate']) for r in values),
            mean_se=statistics.mean(float(r['se']) for r in values),
            coverage=statistics.mean(float(r['ci_lower']) <= float(r['truth']) <= float(r['ci_upper']) for r in values))
        if any(abs(expected[field]-float(row[field])) > 1e-12 for field in expected):
            raise ValueError('Summary does not reproduce its inputs')
    for task in read_rows(study/'manifest.csv'):
        metadata = json.loads((study/'tasks'/task['task_id']/'metadata.json').read_text())
        recorded = metadata['task'][0] if isinstance(metadata['task'],list) else metadata['task']
        repeat = int(task['sim_id'])
        if (recorded['config'] != 'O_moderate_mix' or int(recorded['sim_id']) != repeat or
                int(recorded['population_seed']) != 929000+repeat or int(recorded['fit_seed']) != 939000+repeat):
            raise ValueError('Supplemental seed identity changed')
    mapping = {r['method']:r for r in summary}
    r, t = mapping['RoCE-r3'], mapping['Target-only']
    final_output = output
    output = output.with_name(output.name+'.pending.'+uuid4().hex)
    shutil.copytree(str(previous), str(output), ignore=shutil.ignore_patterns(
        'CHECKS_PASSED','artifact_checksums.json','VISUAL_REVIEW_PASSED','SEED_AUDIT_PASSED'))
    copies = {
        review/'summary.csv':'tables/or_branch_mc200.csv',
        review/'paired_mse.csv':'tables/or_branch_paired_mse.csv',
        audit/'target_law_parity.csv':'tables/or_branch_target_parity.csv',
        study/'reviews/population_projection_v1/population_limits.csv':'tables/or_branch_population_limits.csv',
        base/'model_validation_or_branch_template_v1/population_model_checks.csv':'tables/or_branch_model_checks.csv',
        figures/'rhc_model_mc200.pdf':'figures/or_branch_validation.pdf',
        figures/'figure_preview.png':'figures/or_branch_validation.png',
        figures/'caption.txt':'figures/or_branch_caption.txt'}
    for source, relative_destination in copies.items():
        shutil.copy2(str(source), str(output/relative_destination))
    text = (output/'README_zh.md').read_text()
    text = text.replace('源模型的训练与最终估计使用同一截断规则，并完整重拟合模型和聚合权重。不是只改变最终SE。',
        '源模型的训练与最终估计使用同一截断规则，模型、聚合权重和标准误均据此重新计算。')
    text = text.replace('没有把这项不利结果隐藏或当作已验证的方法。',
        '该探索结果保留在tables中，没有纳入推荐设置。')
    text = text.replace('随后720个新增repeat全部完成，无记录失败。',
        '随后720个新增repeat全部完成，没有未恢复的科学计算失败。')
    text += '''

另完成一个专门核对校准条件的OR正确分支，200 repeats：源协变量分布取90%% target分布加10%%原source经验分布；target分布、真实结局模型、propensity函数、样本量和原有实验均不变。真joint weights约为0.071–9.321，严格处于M=3范围内；log-weight在线性特征基上的投影残差证明权重工作模型错设，OR模型则正确。

该分支的12个初始／校准人口矩系统已直接求解，最大梯度残差约1.36e-8，Hessian最小特征值大于0.00417，截断边界余量至少0.997。它补充了OR正确分支的关键模型与正则性核对；强偏移场景继续作为边界压力测试保留。

M=3的bias为%.3fpp、RMSE为%.3fpp，calibrated Target-only RMSE为%.3fpp，RMSE下降%.1f%%；coverage为%.1f%%（Monte Carlo95%% Wilson区间[%.1f%%,%.1f%%]）。所有200个target-only和oracle结果与原OR场景逐次匹配，差异仅来自source人群分布。

[额外OR分支验证图](figures/or_branch_validation.pdf)及其完整数值位于tables。累计完成五个固定模型各200 repeats；第五个模型不替换、不合并稀释前四个模型，也不改变真实RHC主表。实际临床偏差仍不可由这些实验识别。
''' % (100*float(r['bias']),100*float(r['rmse']),100*float(t['rmse']),
        100*(1-float(r['rmse'])/float(t['rmse'])),100*float(r['coverage']),
        100*float(r['coverage_lower']),100*float(r['coverage_upper']))
    (output/'README_zh.md').write_text(text)
    configuration = json.loads((output/'recommended_configuration.json').read_text())
    configuration['or_branch_validation'] = str(study)
    (output/'recommended_configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    checks = json.loads((output/'audit.json').read_text())
    checks.update(additional_or_branch_repeats=200,total_scenario_repeats=1000,
        additional_or_branch_target_parity=True,additional_or_branch_population_limits_checked=True)
    (output/'audit.json').write_text(json.dumps(checks,indent=2)+'\n')
    provenance = json.loads((output/'provenance.json').read_text())
    provenance.update({str(source):hashlib.sha256(source.read_bytes()).hexdigest() for source in copies})
    (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (output/'CHECKS_PASSED').write_text('Clinical analysis andfive prescribed200-repeat model panels audited;original results andassumption boundaries retained.\n')
    (output/'VISUAL_REVIEW_PASSED').write_text('Original clinical/MC200 figures andsupplemental OR-branch figure were visually reviewed.\n')
    (output/'SEED_AUDIT_PASSED').write_text('The original800 andsupplemental200 scenario/repeat outputs match their prescribed seed identities.\n')
    paths = [p for p in output.rglob('*') if p.is_file() and p.name != 'artifact_checksums.json']
    (output/'artifact_checksums.json').write_text(json.dumps({str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2)+'\n')
    output.rename(final_output)
    print(final_output/'README_zh.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--base',type=Path,default=Path('/scratch.global/zhan9381/FACE-HD/real_data/rhc'))
    args = parser.parse_args()
    append(args.base,args.output)

#!/usr/bin/env python3
"""Compare completed target-specific RHC studies without treating contrasts as bias."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def review(root, output):
    if not str(output.resolve()).startswith('/scratch.global/zhan9381/FACE-HD/real_data/rhc/'):
        raise ValueError('Keep the review in the RHC scratch hierarchy')
    registry = json.loads((root / 'studies.json').read_text())
    assert len(registry) == 4
    assert {study['target_site'] for study in registry} == {'Private', 'Medicare', 'Private & Medicare', 'Medicaid'}
    paired, candidates, baselines = [], [], []
    for study in registry:
        method_root = Path(study['methods'])
        baseline_root = Path(study['baselines'])
        weight_root = Path(study['weight_review'])
        assert (weight_root / 'CHECKS_PASSED').is_file()
        assert (Path(study['paired_review']) / 'CHECKS_PASSED').is_file()
        variance = read_rows(weight_root / 'site_variance.csv')
        weights = read_rows(weight_root / 'source_weight_summary.csv')
        balance = read_rows(weight_root / 'weighted_balance.csv')
        profiles = ['RHC_{}_fold{}'.format(rule, seed)
                    for seed in (0, 101, 202, 303)
                    for rule in ('min', 'source_all_1se', 'global_1se')]
        method_tasks = [row for row in read_rows(method_root / 'manifest.csv') if row['config'] in profiles]
        assert len(method_tasks) == 12
        reference = None
        expected_data_hash = None
        for task in method_tasks:
            directory = method_root / 'tasks' / task['task_id']
            assert (directory / 'COMPLETE').is_file()
            metadata = json.loads((directory / 'metadata.json').read_text())
            assert metadata['target_site'] == study['target_site']
            data_hash = (directory / 'data_sha256.txt').read_text()
            if expected_data_hash is None:
                expected_data_hash = data_hash
            assert expected_data_hash == data_hash
            methods = {row['method']: row for row in read_rows(directory / 'methods.csv')}
            fitted, anchor, ordinary = [methods[name] for name in ('RoCE', 'Calibrated target-only', 'Target-only')]
            if task['config'] == 'RHC_min_fold0':
                reference = ordinary
            profile = task['config']
            site_variance = [row for row in variance if row['profile'] == profile]
            source_weights = [row for row in weights if row['profile'] == profile]
            balance_rows = [row for row in balance if row['profile'] == profile]
            defined = [row for row in balance_rows if row['weighted_minus_target_smd'] != 'NA']
            undefined = [row for row in balance_rows if row['weighted_minus_target_smd'] == 'NA']
            assert all(float(row['target_sd']) == 0 for row in undefined)
            gaps = [abs(float(row['weighted_minus_target_smd'])) for row in defined]
            assert all(math.isfinite(gap) for gap in gaps)
            assert len(site_variance) == 4 and len(source_weights) == 6 and gaps
            assert abs(sum(float(row['variance']) for row in site_variance) - float(fitted['se'])**2) < 1e-12
            paired.append(dict(target_site=study['target_site'], variant=profile[4:].rsplit('_fold', 1)[0],
                fold_seed=int(task['fold_seed']), estimate_pp=100*float(fitted['estimate']),
                se_pp=100*float(fitted['se']), ci_lower_pp=100*float(fitted['ci_lower']),
                ci_upper_pp=100*float(fitted['ci_upper']), target_only_pp=100*float(ordinary['estimate']),
                calibrated_target_pp=100*float(anchor['estimate']),
                gap_to_calibrated_target_pp=100*(float(fitted['estimate'])-float(anchor['estimate'])),
                se_ratio_to_calibrated_target=float(fitted['se'])/float(anchor['se']),
                largest_record_variance_share=max(float(row['top1_total_variance_share']) for row in site_variance),
                min_source_ess_fraction=min(float(row['effective_n'])/float(row['n']) for row in source_weights),
                mean_abs_smd=statistics.mean(gaps), max_abs_smd=max(gaps),
                undefined_smd_count=len(undefined),
                max_constant_feature_mean_gap=max([abs(float(row['weighted_minus_target_mean'])) for row in undefined] or [0.])))
            for source in read_rows(directory / 'sources.csv'):
                label = metadata['source_sites'][int(source['source'][1:])-1]
                candidates.append(dict(target_site=study['target_site'], profile=profile, source_site=label,
                    estimate_pp=100*float(source['estimate']),
                    gap_to_calibrated_target_pp=100*(float(source['estimate'])-float(anchor['estimate'])),
                    weight_mu1=float(source['weight_mu1']), weight_mu0=float(source['weight_mu0']),
                    mean_wald_mu1=float(source['mean_wald_statistic_mu1']),
                    mean_wald_mu0=float(source['mean_wald_statistic_mu0'])))
        assert reference is not None
        baseline_tasks = [row for row in read_rows(baseline_root / 'manifest.csv') if row['role'] == 'baseline']
        assert len(baseline_tasks) == 4
        for task in baseline_tasks:
            directory = baseline_root / 'tasks' / task['task_id']
            assert (directory / 'COMPLETE').is_file()
            assert (directory / 'data_sha256.txt').read_text() == expected_data_hash
            methods = {row['method']: row for row in read_rows(directory / 'methods.csv')}
            for field in ('estimate', 'se'):
                assert abs(float(methods['Target-only'][field])-float(reference[field])) < 1e-10
            fit = methods[task['baseline']]
            baselines.append(dict(target_site=study['target_site'], method=task['baseline'],
                estimate_pp=100*float(fit['estimate']), se_pp=100*float(fit['se']),
                ci_lower_pp=100*float(fit['ci_lower']), ci_upper_pp=100*float(fit['ci_upper'])))
    output.mkdir(parents=True, exist_ok=False)
    write_rows(output / 'paired_results.csv', paired)
    write_rows(output / 'source_candidates.csv', candidates)
    write_rows(output / 'baselines.csv', baselines)
    report = ['# RHC target 轮换诊断', '',
        '同一批5,039名患者，轮换四个保险组为target。每个target对应不同的目标人群和TATE；跨target差异不等于偏差。',
        '所有规定分折均纳入，baselines与各自target的默认分折及普通target参照配对。', '',
        '| Target | 规则 | TATE范围（百分点） | SE范围（百分点） | 相对calibrated target的差值范围（百分点） | 最大病例方差占比 |',
        '|---|---|---:|---:|---:|---:|']
    for study in registry:
        for variant in ('min', 'source_all_1se', 'global_1se'):
            selected = [row for row in paired if row['target_site'] == study['target_site'] and row['variant'] == variant]
            assert len(selected) == 4
            def span(field):
                return '{:.2f}–{:.2f}'.format(min(row[field] for row in selected), max(row[field] for row in selected))
            report.append('| {} | {} | {} | {} | {} | {:.1f}% |'.format(study['target_site'], variant,
                span('estimate_pp'), span('se_pp'), span('gap_to_calibrated_target_pp'),
                100*max(row['largest_record_variance_share'] for row in selected)))
    report += ['', '## Baselines：各target的默认分折', '',
        '| Target | 方法 | 估计（百分点） | SE（百分点） |', '|---|---|---:|---:|']
    report += ['| {} | {} | {:.2f} | {:.2f} |'.format(row['target_site'], row['method'], row['estimate_pp'], row['se_pp'])
               for row in baselines]
    report += ['', '完整数据：[逐分折结果](paired_results.csv)、[source候选及两臂权重](source_candidates.csv)、[baselines](baselines.csv)。', '',
        'target-only也可能有模型错设或未测混杂。候选与target-only的差异是诊断量，不是真实bias估计；不能按哪个target最显著来挑选主分析。',
        '原始min、source-only1se和global1se全部保留。SE、ESS和稳定性改善均不证明消除了混杂或迁移偏差。',
        'target内无变异的特征，其SMD标为NA并另报原始加权均值差；SMD汇总只包含可定义的部分。']
    (output / 'report_zh.md').write_text('\n'.join(report)+'\n')
    (output / 'CHECKS_PASSED').write_text('All four targets,48 method profiles and16 baselines verified; includes reused Private results.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    review(args.root, args.output)

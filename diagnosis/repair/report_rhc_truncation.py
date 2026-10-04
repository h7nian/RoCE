#!/usr/bin/env python3
"""Assemble the audited RHC recommendation and all supporting diagnostics."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def report(base, output):
    base, output = base.resolve(), output.resolve()
    if not str(output).startswith(str(base)+'/'):
        raise ValueError('Keep the report in the RHC scratch hierarchy')
    presentation = base/'presentation_private_min_radius3_v1'
    selected = base/'report_selection_private_min_radius3_v1'
    mc = base/'model_validation_mc200_extension_v1/combined_mc200'
    review = mc/'reviews/complete_v1'
    figures = base/'model_validation_mc200_extension_v1/reviews/figures_v1'
    for directory in (presentation, mc, review, figures):
        if not (directory/'CHECKS_PASSED').is_file():
            raise ValueError('Missing completed audit: '+str(directory))
    methods = read_rows(presentation/'rhc_private_min_radius3_tate.csv')
    by_method = {row['method']:row for row in methods}
    if list(by_method) != ['Target-only','SS','IVW','Federated-DR','Pooled-DR','RoCE']:
        raise ValueError('Unexpected displayed methods or order')
    for row in methods:
        estimate, se = float(row['estimate']), float(row['se'])
        if se <= 0 or any(abs(float(row[field]) - (estimate + sign*1.959963984540054*se)) > 1e-12
                          for field,sign in (('ci_lower',-1),('ci_upper',1))):
            raise ValueError('Invalid reported interval')
    previous = {r['method']:r for r in read_rows(base/'presentation_private_min_calibrated_target_v1/rhc_private_min_tate.csv')}
    for name in set(by_method)-{'RoCE'}:
        for field in ('estimate','se','ci_lower','ci_upper'):
            if float(by_method[name][field]) != float(previous[name][field]):
                raise ValueError('A matched baseline or target anchor changed')
    configuration = json.loads((selected/'configuration.json').read_text())
    if (configuration['source_radius'] != 3 or configuration['target_radius'] != 5 or
            configuration['report_target_site'] != 'Private' or configuration['report_nuisance_lambda_rule'] != 'min'):
        raise ValueError('Unexpected selected analysis')
    summary = read_rows(review/'summary.csv')
    replicates = read_rows(review/'replicates.csv')
    groups = collections.defaultdict(list)
    for row in replicates:
        groups[(row['scenario'],row['label'])].append(row)
    if len(groups) != 20 or len(replicates) != 4000:
        raise ValueError('The full four-by-200 panel is required')
    for row in summary:
        group = groups[(row['scenario'],row['method'])]
        if sorted(int(r['repeat_id']) for r in group) != list(range(1,201)):
            raise ValueError('Missing, repeated or replaced simulation seed')
        errors = [float(r['estimate'])-float(r['truth']) for r in group]
        values = dict(bias=statistics.mean(errors), rmse=math.sqrt(statistics.mean(e*e for e in errors)),
            empirical_sd=statistics.stdev(float(r['estimate']) for r in group),
            mean_se=statistics.mean(float(r['se']) for r in group),
            coverage=statistics.mean(r['covered']=='TRUE' for r in group))
        if any(abs(values[field]-float(row[field])) > 1e-12 for field in values):
            raise ValueError('Monte Carlo summary does not reproduce saved replicates')
    output.mkdir(parents=True, exist_ok=False)
    for name in ('figures','tables'):
        (output/name).mkdir()
    copies = {
        presentation/'rhc_private_min_radius3_tate.csv':'tables/real_data_methods.csv',
        presentation/'rhc_private_min_radius3_sites.csv':'tables/real_data_sites_and_weights.csv',
        presentation/'rhc_private_min_radius3_arms.csv':'tables/real_data_arms.csv',
        presentation/'rhc_private_min_radius3_tate.pdf':'figures/rhc_tate.pdf',
        presentation/'figure_preview.png':'figures/rhc_tate.png',
        figures/'rhc_model_mc200.pdf':'figures/model_validation_mc200.pdf',
        figures/'figure_preview.png':'figures/model_validation_mc200.png',
        figures/'caption.txt':'figures/model_validation_caption.txt',
        review/'summary.csv':'tables/model_validation_mc200.csv',
        review/'paired_mse.csv':'tables/paired_mse_differences.csv',
        base/'source_truncation_confirmation_v1/reviews/report_v1/summary.csv':'tables/fold_stability.csv',
        base/'case_refit_v1/reviews/complete_v1/comparisons.csv':'tables/case_refits.csv',
        base/'reviews/truncation_assessment_20260928_v1/default_partition.csv':'tables/truncation_sensitivity.csv',
        base/'model_validation_mc200_extension_v1/reviews/population_overlap_v1/infeasibility_certificates.csv':
            'tables/population_moment_boundaries.csv',
        base/'mass_stabilization_checks_v1/summary.csv':'tables/exploratory_mass_stabilization.csv'}
    for source, destination in copies.items():
        shutil.copy2(str(source),str(output/destination))
    (output/'recommended_configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    precision_gain = 100*(1-float(by_method['RoCE']['se'])/float(by_method['Target-only']['se']))
    text = '''推荐保留 Private target、原始 min 规则和全部5039位患者，采用 source 截断半径 M=3；target 半径仍为5，aggregation cutoff仍为2。源模型的训练与最终估计使用同一截断规则，并完整重拟合模型和聚合权重。不是只改变最终SE。

当前真实数据结果（百分点）：

| 方法 | TATE | SE | 95% CI |
|---|---:|---:|---:|
'''
    for row in methods:
        text += '| {} | {:.2f} | {:.2f} | [{:.2f}, {:.2f}] |\n'.format(row['method'],
            *[100*float(row[k]) for k in ('estimate','se','ci_lower','ci_upper')])
    text += '''
RoCE的SE比calibrated Target-only低%.1f%%，CI仍包含0。与原M=5的7.90pp、SE4.51pp相比，极端source病例的影响明显减弱。四个baselines和target参照使用相同患者与原定分折，数值未变；图中只显示一个Target-only，代表calibrated anchor。

[RHC对比图](figures/rhc_tate.pdf)沿用Beamer的配色、网格、字体设置和方法顺序。

20组预先指定的新分折中，M=3的TATE范围为2.69–4.76pp，SE为2.11–2.37pp；19/20组SE低于calibrated target-only，中位SE比为0.896。它们是同一队列的分折敏感性检查，不是20份独立临床样本。

18个完整病例重拟合检查全部完成：9个全样本对照复现原结果，9个诊断性删除重拟合全部模型和eta。原分折中，Medicaid病例423删除后，M=5的估计变化为-4.57pp，M=3为-0.55pp。另一分折的对应变化是-6.23与-1.16pp。Target病例303的变化在M=3下也可能更大；完整9个扰动均保留。正式分析不删除这些病例。

四个固定辅助人口模型已各完成200 repeats，三个半径都在相同模拟数据上配对。真实总体TATE由已知人口模型精确求和，不随repeat改变。原20个repeats保留，新增21–200；没有只汇总成功的重复。

| 辅助场景 | M=3 bias (pp) | M=3 RMSE (pp) | Target-only RMSE (pp) | RMSE下降 | M=3 coverage |
|---|---:|---:|---:|---:|---:|
''' % precision_gain
    names = {'O_case_mix':'结局模型正确／较强人群差异','W_overlap':'权重模型正确／较好重叠',
             'W_tail':'权重模型正确／较重尾部','W_departure':'一个source治疗组偏离'}
    lookup = {(r['scenario'],r['method']):r for r in summary}
    for scenario, label in names.items():
        r, t = lookup[(scenario,'RoCE-r3')], lookup[(scenario,'Target-only')]
        text += '| {} | {:.3f} | {:.3f} | {:.3f} | {:.1f}% | {:.1f}% |\n'.format(label,
            100*float(r['bias']),100*float(r['rmse']),100*float(t['rmse']),
            100*(1-float(r['rmse'])/float(t['rmse'])),100*float(r['coverage']))
    text += '''
这些是model-based辅助实验，不是对真实RHC偏差的估计。每组200次的coverage仍有Monte Carlo误差；完整Wilson区间和配对MSE差异见tables。三个W场景共用目标人口模型与配对随机数，不能把它们视为额外独立的target样本。

M=2在强人群差异场景的RMSE略小，但偏差为0.492pp；M=3为0.171pp。真实RHC的平均绝对SMD也从M=2的约0.173改善到M=3的约0.129。因此推荐M=3作为偏差、稳定性与精度的折中，不根据CI是否显著挑选半径。

[辅助验证图](figures/model_validation_mc200.pdf)保留全部半径。误差条是Monte Carlo不确定性：bias使用正态近似，RMSE使用delta近似，coverage使用Wilson区间。

有一个明确的理论边界必须保留：强人群差异的辅助模型中，M=3下Medicaid治疗组的收入参照组校准矩目标值约0.03157，任何受此上界约束的权重最多只能提供0.01556。因此不能仅凭OR模型正确，就声称该场景满足全部校准推断假设。M=2的不可行约束更多；M=5未触发这些逐坐标必要条件，但这也不是全部假设成立的证明。这个人口模型证书不等于已经识别了真实RHC的偏差。

辅助实验中的总体真值、设计概率和oracle方差经过独立核对。额外2000次oracle-only抽样与精确人口方差一致；这些抽样没有混入完整拟合的200-repeat比较。

数值可靠性方面，验证暴露并修复了初始／校准权重CV网格上界不足的问题。只有整条原网格失败时，才在相同CV折追加更大的lambda；原候选、min规则、优化器和收敛标准保留。20项回归测试、313个断言通过，两份真实失败输入恢复。原成功路径的12个辅助拟合和真实RHC M=3/M=5完整复现通过，点估计及SE差值为0或在1e-10容差内。

原80个pilot任务曾有8个校准CV失败，均在原种子下完整重跑；另4个成功repeat作一致性对照。随后720个新增repeat全部完成，无记录失败。最终四场景各200，所有旧失败与旧结果保留。

另试过保留M=5、按source权重总量收缩聚合增量的研究方案。原分折得到4.02pp、first-order SE2.53pp，未优于target-only的2.46pp；它也需要额外的推断验证，因此不作为本次推荐。没有把这项不利结果隐藏或当作已验证的方法。

本次交付是可复核的分析与图表结果包，没有向Overleaf推送。原M=5、其他截断半径、source1se和其他target的旧结果都保留；本报告的真实数据主表仅采用Private/min。
'''
    (output/'README_zh.md').write_text(text)
    provenance = {str(source):hashlib.sha256(source.read_bytes()).hexdigest() for source in copies}
    provenance[str(selected/'configuration.json')] = hashlib.sha256((selected/'configuration.json').read_bytes()).hexdigest()
    (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (output/'audit.json').write_text(json.dumps(dict(
        cohort_and_partition_unchanged=True, nuisance_min_rule_preserved=True,
        source_training_and_evaluation_radius=3,target_radius=5,aggregation_cutoff=2,
        only_calibrated_target_displayed=True,matched_baselines_unchanged=True,
        all_800_prescribed_simulation_repeats_present=True,metrics_recomputed_from_4000_rows=True,
        clinical_truth_known=False,population_calibration_boundary_disclosed=True,
        clinical_interval_includes_zero=True,old_results_preserved=True,
        manuscript_pushed=False),indent=2)+'\n')
    (output/'CHECKS_PASSED').write_text('Selected clinical results,unchanged comparators andcomplete MC200 metrics checked; limitations explicitly retained.\n')
    print(output/'README_zh.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--base',type=Path,default=Path('/scratch.global/zhan9381/FACE-HD/real_data/rhc'))
    args = parser.parse_args()
    report(args.base,args.output)

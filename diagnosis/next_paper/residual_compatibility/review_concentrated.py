"""Report all paired calibration ablations and the equally improved target check."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = ['#777777', '#0072B2', '#D55E00', '#009E73']


def read_rows(path):
    return list(csv.DictReader(path.open()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    args = parser.parse_args()
    root = Path(args.root)
    report = root / 'report'
    report.mkdir(exist_ok=True)
    binary = read_rows(root / 'binary/metrics.csv')
    fitted = read_rows(root / 'fitted/metrics.csv')
    target = read_rows(root / 'paired_target_audit/metrics.csv')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'axes.grid': True, 'grid.alpha': .2, 'grid.linewidth': .5})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for axis, floor in zip(axes, ['0.0', '0.3']):
        for method, label, color in zip(['hybrid_legacy', 'fourier_legacy', 'hybrid_direct', 'dispersion_direct'],
                ['Previous hybrid', 'Fourier only', 'Direct hybrid', 'Direct dispersion'], COLORS):
            rows = sorted([r for r in binary if r['pattern'] == 'weak_boundary' and r['shared_scale'] == floor and r['method'] == method], key=lambda r: int(r['source_count']))
            axis.plot([int(r['source_count']) for r in rows], [float(r['length_ratio_target']) for r in rows], 'o-', color=color, label=label)
        axis.axhline(1, color='black', ls=':', lw=1)
        axis.set_xscale('log', base=2)
        axis.set_xticks([32, 512, 2048])
        axis.set_xticklabels(['32', '512', '2048'])
        axis.set_xlabel('Number of sources K')
        axis.set_title('Shared target heterogeneity = ' + floor)
    axes[0].set_ylabel('Mean length / ordinary target length')
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(str(report / 'bounded_comparison.pdf'))
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for method, label, color in zip(['hybrid_legacy', 'hybrid_aggregate_budget', 'dispersion_aggregate_budget'],
            ['Previous hybrid', 'Aggregate-budget hybrid', 'Aggregate-budget dispersion'], [COLORS[0], COLORS[2], COLORS[3]]):
        rows = sorted([r for r in fitted if r['pattern'] == 'weak_boundary' and r['shared_scale'] == '0.3' and r['method'] == method], key=lambda r: int(r['source_count']))
        ratios = []
        for row in rows:
            denominator = next(float(t['mean_length']) for t in target if t['cell_id'] == row['cell_id'] and t['method'] == 'full_target_direct_budget')
            ratios.append(float(row['mean_length']) / denominator)
        axes[0].plot([int(r['source_count']) for r in rows], ratios, 'o-', label=label, color=color)
    axes[0].axhline(1, color='black', ls=':', lw=1)
    axes[0].set_ylabel('Length / improved full target length')
    axes[0].legend(frameon=False, fontsize=8)
    for key, label, color in [('old_source_budget', 'Previous source budget', COLORS[0]),
                              ('direct_source_budget', 'Aggregate source budget', COLORS[3]),
                              ('oracle_worst_subset_mean', 'Oracle worst subset mean', COLORS[1])]:
        values = []
        for cell in [4, 8, 12]:
            diagnostic = json.loads((root / 'fitted/diagnostics' / ('cell_%03d.json' % cell)).read_text())
            values.append(np.mean(diagnostic[key]))
        axes[1].plot([8, 32, 128], values, 'o-', label=label, color=color)
    axes[1].set_ylabel('Mean nuisance allowance')
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xscale('log', base=2)
        axis.set_xticks([8, 32, 128])
        axis.set_xticklabels(['8', '32', '128'])
        axis.set_xlabel('Number of sources K')
    fig.tight_layout()
    fig.savefig(str(report / 'fitted_comparison.pdf'))
    plt.close(fig)
    # Paired SEs preserve the covariance between procedures on the same data.
    comparisons = []
    for profile, rows in [('binary', binary), ('fitted', fitted)]:
        for cell in sorted(set(int(r['cell_id']) for r in rows)):
            packet = np.load(str(root / profile / 'cells' / ('cell_%03d.npz' % cell)))
            names = list(packet['methods'])
            lengths = np.diff(packet['interval'], axis=2)[:, :, 0]
            reference = lengths[:, names.index('hybrid_legacy')]
            for index, name in enumerate(names):
                difference = lengths[:, index] - reference
                comparisons.append(dict(profile=profile, cell_id=cell, method=name,
                    mean_length_difference=float(difference.mean()), mcse=float(difference.std(ddof=1)/np.sqrt(len(difference)))))
    with (report / 'paired_differences.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    summary = dict(binary_settings=18, binary_repeats=1000, fitted_settings=12, fitted_repeats=300,
        total_setting_repetitions=21600, all_results_retained=True,
        minimum_observed_coverage=min(float(r['coverage']) for r in binary + fitted),
        aggregate_certificate_failure=max(float(r['certificate_failure']) for r in target),
        note='The original fitted diagnostic flagged machine-roundoff in the exactly zero d drift. The paired_target_audit is authoritative and uses tolerance 1e-12; intervals are unchanged.')
    (report / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (report / 'REPORT_zh.txt').write_text('''RoCE-K v16：整体误差校准与配对消融，2026-10-02

主线仍是 tau=d+r1-r0：target 约束共同预测与两臂残差，source 检验候选残差中心。
不要求先筛出每个弱偏离来源。新增的是整体误差控制，不是新的聚合权重。

完成：18 个有界 oracle-score 设置，各1000 repeats；12 个重新拟合四-strata设置，
各300 repeats，共21600个setting-specific重复。每站总样本量1000。相同数据、相同总
错误概率预算配对比较 Fourier-only、dispersion-only、hybrid；新旧结果都保留。
新的随机数种子与 v15 不同，不能将两版独立样本差值当作配对差值。

1. 推导并实现患者级零对角二次统计量的指数型上界，直接处理离散度，
   无需逐source方差证书。来源平均值用允许独立非同分布患者的 empirical Bernstein。
   完成4096原子非高斯精确枚举、矩阵谱检查、独立LP求解器核对及相关测试。
2. 四strata中利用source独立pilot的校准投影、target概率单纯形约束，以及outcome误差
   box的角点，保护每个未知g-subset的有符号平均余项和平方余项。不是直接把最大值改平均。
3. 区分已知方差修正D-circ和随机方差修正D-hat；只有前者的Gaussian同方差比较具有
   精确非中心卡方分布。历史有效的代码路径使用了正确统计量；修复稿件记号并加API保护。

实证结果：
- binary、K=2048、weak bias、shared heterogeneity=0.3：旧hybrid长度为普通target
  normal区间的1.736倍；新hybrid为1.079倍；单用直接离散度为1.066倍。
- 没有共享target异质性时，同一weak-bias设置的新离散度比例为0.301，旧hybrid为0.823。
- 四strata重新拟合、K=128、weak bias、heterogeneity=0.3：treated/control平均source
  预算从0.285/0.263降至0.111/0.097。真实最坏g-subset平均漂移仍仅约0.0034/0.0032。
- 同一拟合设置，旧hybrid平均长度1.110；新hybrid0.882；新离散度0.872。
  但同步改进预算后的全target两折区间长度0.831。新方法仍未获得净借用优势。
  因此不能只用旧的target证书作为比较，宣称已实现实用效率提升。
- 当前所有组件方法的观察coverage均为1；这说明仍保守，不是覆盖率恰好95%的证据。
  每设置1000次全覆盖时二侧95% Clopper-Pearson下界约0.9963，300次时约0.9878。
- 当前有界分数保护下，Fourier加入hybrid没有显示净长度优势；它额外占用错误预算。
  保留该选项及Gaussian理想比较，不据此宣称Fourier在所有模型中无用。

重要诊断修正：fitted/diagnostics内最初certificate_failure把理论上为0的d漂移与
浮点数(约1e-16)作严格比较，产生伪失败。保留该原记录；paired_target_audit按1e-12
容差重新检查，所有3600重复均无实质证书失败。任何估计、区间和coverage结果均未改动。

边界：这些是有界oracle分数和四strata诚实拆分实验，不是高维完整RoCE的证明。
固定strata的角点算法对K线性，但对strata数指数；不能直接推广成高维实现。
K增长的有限样本覆盖与区间变短是两个结论；共享target噪声、拟合余项和Fourier
非高斯/方差校正仍可限制效率。所有失败与fallback均保留，没有更改DGP以隐藏问题。

下一步：优先继续压缩整体nuisance余项，对比同样改进的全target方法；再证明可扩展的
校准特征版本。直接离散度是当前有界分数研究中更有用的简洁候选；Fourier作为
预先指定的增强/消融保留。现有生产RoCE/ENAR算法与旧实验均未替换。
''')


if __name__ == '__main__':
    main()

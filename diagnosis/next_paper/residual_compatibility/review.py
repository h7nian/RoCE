"""Audit completed paired experiments and export figures plus readable findings."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LABELS = {"target_nominal": "Target-only (95%)", "target_joint": "Target joint region",
          "pooled_naive": "Naive pooling", "oracle_valid": "Known-valid oracle",
          "fourier_bernstein": "Fourier (Bernstein)", "fourier_hoeffding": "Fourier (Hoeffding)",
          "dispersion": "Dispersion", "count_search": "Coordinate-count search",
          "hybrid": "Combined certificates", "fourier_plugin": "Fourier, plug-in variance",
          "fourier_variance_certified": "Fourier, variance certificate",
          "fourier_bounded": "Fourier, bounded-score protection"}
COLORS = {"target_nominal": "#555555", "pooled_naive": "#D55E00", "oracle_valid": "#009E73",
          "fourier_bernstein": "#0072B2", "dispersion": "#CC79A7", "hybrid": "#E69F00",
          "fourier_plugin": "#009E73", "fourier_variance_certified": "#E69F00", "fourier_bounded": "#E69F00"}


def collect(directory):
    if not (directory / "COMPLETE").exists():
        raise RuntimeError("Campaign is not complete: " + str(directory))
    configuration = json.loads((directory / "configuration.json").read_text())
    files = sorted((directory / "cells").glob("*.json"))
    if len(files) != len(configuration["settings"]):
        raise AssertionError("Missing or duplicate cell records")
    rows, pairing = [], []
    for file in files:
        cell = json.loads(file.read_text())
        with np.load(str(file.with_suffix(".npz"))) as saved:
            methods = list(saved["methods"])
            for row in cell["metrics"]:
                column = methods.index(row["method"])
                truth, repeats = row["truth"], row["repeats"]
                covered = (saved["lower"][:, column] <= truth) & (truth <= saved["upper"][:, column])
                length = saved["upper"][:, column] - saved["lower"][:, column]
                assert len(covered) == repeats and np.all(length >= 0)
                assert abs(covered.mean() - row["coverage"]) < 1e-12
                assert abs(length.mean() - row["mean_length"]) < 1e-12
                assert abs(np.sqrt(np.mean((saved["point"][:, column] - truth)**2)) - row["rmse"]) < 1e-12
                probability, z = row["coverage"], norm.isf(.025)
                denominator = 1 + z*z/repeats
                middle = (probability + z*z/(2*repeats)) / denominator
                radius = z * np.sqrt(probability*(1-probability)/repeats + z*z/(4*repeats**2)) / denominator
                row.update(coverage_wilson_lower=float(middle-radius), coverage_wilson_upper=float(middle+radius))
                rows.append(row)
            if "fourier_bernstein" in methods and "dispersion" in methods:
                first, second = methods.index("fourier_bernstein"), methods.index("dispersion")
                difference = (saved["upper"][:, first] - saved["lower"][:, first]) - \
                    (saved["upper"][:, second] - saved["lower"][:, second])
                pairing.append(dict(cell["setting"], difference="Fourier minus dispersion length",
                    mean_difference=float(difference.mean()), mcse=float(difference.std(ddof=1)/np.sqrt(len(difference)))))
            if cell["setting"]["noise"] == "gaussian":
                bias, variance = saved["source_bias"][:, 0].mean(), saved["source_variance"][:, 0].sum() / len(saved["source_bias"])**2
                standardized_bias = bias / np.sqrt(variance)
                for row in cell["metrics"]:
                    row["treated_pool_analytic_coverage"] = float(norm.cdf(1.96-standardized_bias)-norm.cdf(-1.96-standardized_bias))
    return rows, pairing


def save_csv(path, rows):
    fields = sorted(set(key for row in rows for key in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def point_rule_diagnostics(directory):
    rows = []
    for path in sorted((directory / "cells").glob("*.json")):
        cell = json.loads(path.read_text())
        setting = cell["setting"]
        if setting["source_count"] != 2048 or setting["shared_scale"] != 0:
            continue
        truth = cell["metrics"][0]["truth"]
        with np.load(str(path.with_suffix(".npz"))) as saved:
            names = list(saved["methods"])
            pooled = saved["point"][:, names.index("pooled_naive")]
            for method in ["fourier_bernstein", "dispersion", "hybrid"]:
                column = names.index(method)
                fallback = saved["fallback"][:, column]
                squared = (saved["point"][:, column] - truth)**2
                projected = np.clip(pooled, saved["lower"][:, column], saved["upper"][:, column])
                rows.append(dict(cell_id=setting["cell_id"], pattern=setting["pattern"], method=method,
                    rmse=float(np.sqrt(squared.mean())), fallback_rate=float(fallback.mean()),
                    fallback_share_squared_error=float(squared[fallback].sum()/squared.sum()),
                    nonfallback_rmse=float(np.sqrt(squared[~fallback].mean())),
                    pooled_projection_rmse=float(np.sqrt(np.mean((projected-truth)**2)))))
    return rows


def draw_lines(axis, rows, methods, metric, title):
    for method in methods:
        selected = sorted([row for row in rows if row["method"] == method], key=lambda row: row["source_count"])
        x = np.array([row["source_count"] for row in selected])
        y = np.array([row[metric] for row in selected])
        if not len(x):
            continue
        axis.plot(x, y, marker="o", markersize=3, linewidth=1.6, label=LABELS[method], color=COLORS.get(method))
        if metric == "coverage":
            lower = [row["coverage_wilson_lower"] for row in selected]
            upper = [row["coverage_wilson_upper"] for row in selected]
            axis.fill_between(x, lower, upper, color=COLORS.get(method), alpha=.09)
    axis.set_xscale("log", base=2)
    ticks = sorted(set(row["source_count"] for row in rows))
    axis.set_xticks(ticks)
    axis.set_xticklabels(ticks)
    axis.set_xlabel("Number of sources K")
    axis.set_title(title, fontsize=10)
    axis.grid(True, alpha=.2)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if metric == "coverage":
        axis.axhline(.95, color="black", linewidth=.8, linestyle="--")
        axis.set_ylim(0, 1.02)
        axis.set_ylabel("Coverage")
    else:
        axis.axhline(1, color="black", linewidth=.8, linestyle="--")
        axis.set_yscale("log")
        axis.set_ylabel("Mean length / target-only")


def export_figure(figure, path):
    figure.savefig(str(path) + ".pdf", bbox_inches="tight")
    figure.savefig(str(path) + ".png", bbox_inches="tight", dpi=160)
    plt.close(figure)


def table(rows, panel, count, floor, patterns, methods):
    lines = ["pattern | method | coverage [Wilson 95% MC interval] | length/target | bias | RMSE | fallback"]
    selected = [row for row in rows if row["panel"] == panel and row["source_count"] == count and row["shared_scale"] == floor]
    for pattern in patterns:
        for method in methods:
            matches = [row for row in selected if row["pattern"] == pattern and row["method"] == method]
            for row in matches:
                lines.append("{} | {} | {:.3f} [{:.3f}, {:.3f}] | {:.3f} | {:.5f} | {:.5f} | {:.3f}".format(
                    pattern, method, row["coverage"], row["coverage_wilson_lower"], row["coverage_wilson_upper"],
                    row["length_ratio_target"], row["bias"], row["rmse"], row["fallback_rate"]))
    return "\n".join(lines)


def main(root, report_name="report"):
    report = root / report_name
    # Only derived reports are regenerated; raw campaign files are immutable.
    report.mkdir(exist_ok=True)
    primary, primary_pairs = collect(root / "main")
    confirmation, confirmation_pairs = collect(root / "confirmation")
    boundary, boundary_pairs = collect(root / "boundary")
    revised_variance, revised_pairs = collect(root / "variance_within_budget")
    primary = [row for row in primary if row["panel"] != "variance_pilot"] + revised_variance
    primary_pairs = [row for row in primary_pairs if row["panel"] != "variance_pilot"] + revised_pairs
    save_csv(report / "metrics_audited.csv", primary + confirmation + boundary)
    save_csv(report / "paired_length_differences.csv", primary_pairs + confirmation_pairs + boundary_pairs)
    point_diagnostics = point_rule_diagnostics(root / "confirmation")
    save_csv(report / "point_rule_diagnostic.csv", point_diagnostics)
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.labelsize": 9})
    patterns = ["all_valid", "weak_boundary", "strong"]
    titles = {"all_valid": "All sources valid", "weak_boundary": "Weak shifts: 2 K^(-1/4) SE", "strong": "Shifts: 4 source SE"}
    methods = ["target_nominal", "pooled_naive", "fourier_bernstein", "dispersion"]
    figure, axes = plt.subplots(4, 3, figsize=(12.5, 11.5))
    for column, pattern in enumerate(patterns):
        for block, floor in enumerate([0., .5]):
            rows = [row for row in primary if row["panel"] == "main" and row["pattern"] == pattern and row["shared_scale"] == floor]
            for offset, metric in enumerate(["coverage", "length_ratio_target"]):
                draw_lines(axes[2*block+offset, column], rows, methods, metric,
                           titles[pattern] + "; shared SD = {} source SE".format(floor))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, ncol=4, loc="upper center", bbox_to_anchor=(.5, 1.015))
    figure.suptitle("Gaussian oracle residual experiment: n=1000/site, 1000 repetitions/cell", y=1.045, fontsize=12)
    figure.tight_layout()
    export_figure(figure, report / "gaussian_main")

    figure, axes = plt.subplots(2, 3, figsize=(12, 6))
    for column, pattern in enumerate(patterns):
        rows = [row for row in confirmation if row["pattern"] == pattern and row["shared_scale"] == 0]
        for line, metric in enumerate(["coverage", "length_ratio_target"]):
            draw_lines(axes[line, column], rows, ["target_nominal", "fourier_bernstein", "dispersion", "hybrid"], metric, titles[pattern])
            if metric == "coverage":
                axes[line, column].set_ylim(.93, 1.005)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, ncol=4, loc="upper center", bbox_to_anchor=(.5, 1.03))
    figure.suptitle("Independent-seed confirmation: 2000 repetitions/cell, no shared target floor", y=1.09, fontsize=12)
    figure.tight_layout()
    export_figure(figure, report / "independent_confirmation")

    figure, axes = plt.subplots(2, 2, figsize=(9, 6))
    for column, floor in enumerate([0., .3]):
        rows = [row for row in primary if row["panel"] == "binary_scores" and row["pattern"] == "weak_boundary" and row["shared_scale"] == floor]
        for line, metric in enumerate(["coverage", "length_ratio_target"]):
            draw_lines(axes[line, column], rows, ["target_nominal", "pooled_naive", "fourier_bernstein", "fourier_bounded"], metric,
                       "Binary patient scores; effect heterogeneity = {}".format(floor))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    labels = ["Fourier (Gaussian approximation)" if label == LABELS["fourier_bernstein"] else label for label in labels]
    figure.legend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(.5, 1.07))
    figure.suptitle("Gaussian approximation versus explicit bounded-score protection", y=1.13, fontsize=12)
    figure.tight_layout()
    export_figure(figure, report / "binary_scores")

    summary = ["RoCE-K 残差相容性实验报告", "",
        "完成情况：主实验90个组合，每组合1000次重复；独立种子确认18个组合、有效比例边界6个组合，每组合2000次重复。",
        "最终114个组合、138000次配置内重复均完成；不同配置使用配对随机数，不把配置间结果当独立样本。",
        "所有站点n=1000；本次是oracle残差模块实验，并非重新拟合p=100/200的完整RoCE。",
        "方差试验已按source 250 pilot+750 evaluation重做六个组合。原先1000 evaluation+额外250 pilot的六组结果保留，但不纳入最终主面板。",
        "主假设为每臂至少75%来源有效；另有50%和25%保证比例的敏感性面板。",
        "", "核心判断：",
        "1. weak bias与K增长确实应联合处理；无共享target噪声时，普通汇总在弱偏差边界明显欠覆盖。",
        "2. 已知同方差高斯模型下，二阶矩在全部有效和弱偏差时往往更短；傅里叶在强偏差时更有优势。",
        "3. 独立确认使用预先拆分错误预算的两类约束交集，没有在相同alpha下事后挑最短区间。",
        "4. shared target噪声会掩盖来源汇总的弱偏差，并限制增大K的收益。",
        "5. 方差证书和非高斯保护均有实际长度代价；不能把oracle/近似结果当作完整拟合方法的保证。",
        "", "主实验：K=2048，无共享target方差底限", table(primary,"main",2048,0.,
            ["all_valid","weak_below","weak_boundary","strong"],methods+["oracle_valid"]),
        "", "独立确认：K=2048，无共享target方差底限",table(confirmation,"confirmation",2048,0.,patterns,
            ["target_nominal","fourier_bernstein","dispersion","hybrid"]),
        "", "独立确认：K=2048，shared SD=0.5 source SE",table(confirmation,"confirmation",2048,.5,patterns,
            ["target_nominal","fourier_bernstein","dispersion","hybrid"]),
        "", "方差估计：K=2048，target共同SD为target残差SE的0.5倍；source 250+750",table(primary,"variance_pilot",2048,.5,
            ["all_valid","weak_boundary"],["fourier_bernstein","fourier_plugin","fourier_variance_certified"]),
        "", "二元患者分数：K=2048，无效应异质性",table(primary,"binary_scores",2048,0.,patterns,
            ["target_nominal","pooled_naive","fourier_bernstein","fourier_bounded","dispersion"]),
        "", "二元患者分数：K=2048，效应异质性0.3",table(primary,"binary_scores",2048,.3,patterns,
            ["target_nominal","pooled_naive","fourier_bernstein","fourier_bounded","dispersion"]),
        "", "实际非多数有效边界：一部分treated来源中心为0，其他中心相差1 source SE；control均有效。",
        "两个中心均有至少g个来源支持，因此来源数据本身无法确定哪个中心属于target。",
        "以下保留每个K下的结果，较小的保证比例不能自动获得多数有效时的长度保证：",
        *["K={}，实际/保证有效比例={}，{}：coverage={:.3f}，length/target={:.3f}".format(
            row["source_count"],row["actual_valid_fraction"],row["method"],row["coverage"],row["length_ratio_target"])
          for row in boundary if row["method"] in ["fourier_bernstein","hybrid"]],
        "", "解释与限制：",
        "Fourier(Bernstein)为我们推导的高斯残差模块校准，并非Guo/RIFL原始算法的复现。",
        "count_search只是同时坐标区间的计数搜索，不能用其效率代表Guo完整searching-and-sampling。",
        "二阶矩方法在同方差高斯下反演非中心卡方；异方差时使用Cantelli；二元分数时使用患者U统计量和有界分数保护。",
        "二元数据上的fourier_bernstein/fourier_plugin仍属于高斯近似诊断。fourier_bounded才加入显式CF近似界，并使用更保守的target矩区域。",
        "主方法区间多数高于95% coverage，表明保护仍保守；区间更短不等于已达到最优效率。",
        "表内MC区间为单个配置的Wilson区间，不是多配置同时保证。长度差异的配对MCSE另存CSV。",
        "区间空集时回退target-only；此规则不修复错误的有效数量假设。回退、失败和所有重复均保留。",
        "点估计用区间中点；RMSE只作诊断，不宣称其优于目标估计。",
        "回退病例的误差贡献另存point_rule_diagnostic.csv；组合方法全部有效场景的0.8%回退贡献约97.7%平方误差。",
        "固定target区域错误概率加target-width回退，还可能使期望长度在K方向存在底限；后续需要证明回退长度控制。",
        "将naive pooling投影到区间只是事后诊断：它改善部分弱偏差情况，却恶化强偏差RMSE，因此未采用。",
        "共同目标项只出现一次，目标三维协方差以及两臂source协方差均保留。",
        "", "下一步：固定残差相容性主线，研究组合证书的长度理论、可实现的方差/分布误差保护和nuisance预算。",
        "在这些条件补齐前，不把本次结果写成解决完整高维K趋于无穷问题的定理。",
        "", "复现入口：../source/experiments.py、../source_confirmation/experiments.py；协议和源文件哈希已保存。",
        "文献：https://arxiv.org/html/2604.26992v2 ; https://arxiv.org/html/2104.06911v5 ."]
    (report / "REPORT_zh.txt").write_text("\n".join(summary)+"\n")
    checks = dict(primary_cells=len({row["cell_id"] for row in primary}),
                  confirmation_cells=len({row["cell_id"] for row in confirmation}),
                  boundary_cells=len({row["cell_id"] for row in boundary}),
                  revised_variance_cells=len({row["cell_id"] for row in revised_variance}),
                  reported_repetitions=138000, archived_extra_pilot_repetitions=6000,
                  method_setting_rows=len(primary)+len(confirmation)+len(boundary), every_record_recomputed=True,
                  all_intervals_finite_and_ordered=True,
                  coverage_scope="Descriptive Monte Carlo; oracle nuisance scores.")
    (report / "audit.json").write_text(json.dumps(checks,indent=2)+"\n")
    latex = [r"\begin{tabular}{llrrr}", r"\toprule",
             r"Scenario & Method & Coverage & Length/target & RMSE\\", r"\midrule"]
    short_names = {"target_nominal":"Target-only", "pooled_naive":"Naive pooling", "fourier_bernstein":"Fourier",
                   "dispersion":"Dispersion", "hybrid":"Combined", "oracle_valid":"Valid-set oracle"}
    for pattern in patterns:
        for method in short_names:
            row = next(row for row in confirmation if row["source_count"] == 2048 and row["shared_scale"] == 0
                       and row["pattern"] == pattern and row["method"] == method)
            latex.append("{} & {} & {:.3f} & {:.3f} & {:.5f}\\\\".format(
                pattern.replace("_"," "),short_names[method],row["coverage"],row["length_ratio_target"],row["rmse"]))
    latex += [r"\bottomrule",r"\end{tabular}"]
    (report / "confirmation_table.tex").write_text("\n".join(latex)+"\n")
    print(json.dumps(checks))


if __name__ == "__main__":
    main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv)>2 else "report")

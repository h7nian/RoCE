"""Audit and report the fallback/calibration/nuisance follow-up experiments."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LABELS={"fixed_target":"Original target fallback","fixed_shorter":"Shorter protected fallback",
        "decay_target":"Decreasing target error, original fallback",
        "decay_shorter":"Decreasing target error, shorter fallback",
        "known_variance_moment":"Known variance: moment target",
        "known_variance_bernstein":"Known variance: Bernstein target",
        "known_variance_empirical":"Known variance: empirical Bernstein",
        "estimated_variance_empirical":"Estimated variance: full protection",
        "gaussian_plugin":"Gaussian plug-in (diagnostic)",
        "gaussian_zero_budget":"Gaussian, zero bias budget",
        "gaussian_oracle_budget":"Gaussian, oracle bias budget",
        "bounded_oracle_budget":"Bounded scores, oracle bias budget",
        "fully_feasible":"Feasible simultaneous budgets",
        "feasible_shared_target_budget":"Shared-target budget paid once"}
COLORS=["#777777","#0072B2","#E69F00","#009E73","#CC79A7"]


def audit_campaign(path):
    assert (path/"COMPLETE").exists(),str(path)
    plan=json.loads((path/"configuration.json").read_text())
    files=sorted((path/"cells").glob("*.json"))
    assert len(files)==len(plan["settings"])
    rows=[]
    for file in files:
        cell=json.loads(file.read_text())
        with np.load(str(file.with_suffix(".npz"))) as raw:
            for column,row in enumerate(cell["metrics"]):
                assert raw["methods"][column]==row["method"]
                interval=raw["interval"][:,column]
                assert np.all(np.isfinite(interval)) and np.all(interval[:,1]>=interval[:,0])
                covered=(interval[:,0]<=raw["truth"])&(raw["truth"]<=interval[:,1])
                assert abs(covered.mean()-row["coverage"])<1e-12
                assert abs(np.mean(interval[:,1]-interval[:,0])-row["mean_length"])<1e-12
                for field in ["point_midpoint","point_source"]:
                    key="rmse_midpoint" if field=="point_midpoint" else "rmse_source"
                    assert abs(np.sqrt(np.mean((raw[field][:,column]-raw["truth"])**2))-row[key])<1e-12
                if row["method"] not in ["fixed_target","decay_target"]:
                    assert np.all(interval[:,1]-interval[:,0]<=raw["length_envelope"][:,column]+1e-8)
                z=norm.isf(.025);n=len(covered);p=covered.mean()
                middle=(p+z*z/(2*n))/(1+z*z/n)
                radius=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
                row.update(coverage_lower=float(middle-radius),coverage_upper=float(middle+radius))
                rows.append(row)
    return rows


def save_csv(path,rows):
    with path.open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted(set(key for row in rows for key in row)))
        writer.writeheader();writer.writerows(rows)


def panel(axis,rows,methods,metric,title):
    for color,method in zip(COLORS,methods):
        chosen=sorted([r for r in rows if r["method"]==method],key=lambda r:r["source_count"])
        if not chosen:continue
        axis.plot([r["source_count"] for r in chosen],[r[metric] for r in chosen],"o-",markersize=4,
                  label=LABELS[method],color=color,linewidth=1.7)
    axis.set_xscale("log",base=2)
    ticks=sorted(set(r["source_count"] for r in rows));axis.set_xticks(ticks);axis.set_xticklabels(ticks)
    axis.set_xlabel("Sources K");axis.set_title(title,fontsize=10)
    axis.grid(alpha=.2);axis.spines["top"].set_visible(False);axis.spines["right"].set_visible(False)
    if metric=="coverage":
        axis.set_ylim(.94,1.005);axis.axhline(.95,color="black",ls="--",lw=.8);axis.set_ylabel("Coverage")
    else:
        axis.set_yscale("log")
        axis.set_ylabel("Midpoint RMSE" if metric=="rmse_midpoint" else "Length / target-only normal CI")
        if metric=="length_ratio_target":axis.axhline(1,color="black",ls="--",lw=.8)


def save_figure(figure,path):
    figure.tight_layout()
    figure.savefig(str(path)+".pdf",bbox_inches="tight")
    figure.savefig(str(path)+".png",bbox_inches="tight",dpi=160)
    plt.close(figure)


def table(rows,count,floor,pattern,methods):
    output=["method | coverage | length/normal target | midpoint RMSE | source-point RMSE | fallback"]
    for method in methods:
        row=next(r for r in rows if r["source_count"]==count and r["shared_scale"]==floor
                 and r["pattern"]==pattern and r["method"]==method)
        output.append("{} | {:.4f} | {:.3f} | {:.5f} | {:.5f} | {:.4f}".format(method,row["coverage"],
            row["length_ratio_target"],row["rmse_midpoint"],row["rmse_source"],row["fallback_rate"]))
    return "\n".join(output)


def main(root):
    report=root/"report";report.mkdir(exist_ok=True)
    gaussian=audit_campaign(root/"gaussian")
    binary=audit_campaign(root/"binary")
    fitted=audit_campaign(root/"fitted_shared_budget")
    # This refinement uses the same fitted/evaluation data. Its first five
    # method outputs must reproduce the earlier panel, not create a new sample.
    for file in sorted((root/"fitted"/"cells").glob("*.npz")):
        with np.load(str(file)) as before,np.load(str(root/"fitted_shared_budget"/"cells"/file.name)) as after:
            for key in ["interval","point_midpoint","point_source","fallback"]:
                np.testing.assert_allclose(before[key],after[key][:,:len(before["methods"])],atol=1e-12,rtol=0)
    save_csv(report/"metrics_audited.csv",gaussian+binary+fitted)
    plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans"})
    figure,axes=plt.subplots(2,2,figsize=(10,6))
    for column,floor in enumerate([0.,.5]):
        chosen=[r for r in gaussian if r["pattern"]=="weak_boundary" and r["shared_scale"]==floor]
        for line,metric in enumerate(["rmse_midpoint","length_ratio_target"]):
            panel(axes[line,column],chosen,["fixed_target","fixed_shorter","decay_shorter"],metric,
                  "Weak bias; shared target SD = {} source SE".format(floor))
    handles,labels=axes[0,0].get_legend_handles_labels()
    figure.legend(handles,labels,ncol=1,loc="upper center",bbox_to_anchor=(.5,1.13))
    save_figure(figure,report/"fallback_comparison")
    figure,axes=plt.subplots(2,2,figsize=(10,6))
    methods=["known_variance_moment","known_variance_bernstein","known_variance_empirical","estimated_variance_empirical"]
    for column,floor in enumerate([0.,.3]):
        chosen=[r for r in binary if r["pattern"]=="weak_boundary" and r["shared_scale"]==floor]
        for line,metric in enumerate(["coverage","length_ratio_target"]):
            panel(axes[line,column],chosen,methods,metric,"Binary scores; effect heterogeneity = {}".format(floor))
    handles,labels=axes[0,0].get_legend_handles_labels()
    figure.legend(handles,labels,ncol=2,loc="upper center",bbox_to_anchor=(.5,1.1))
    save_figure(figure,report/"binary_calibration")
    figure,axis=plt.subplots(figsize=(8,4))
    chosen=[r for r in fitted if r["pattern"]=="weak_boundary" and r["shared_scale"]==.3]
    panel(axis,chosen,["gaussian_oracle_budget","bounded_oracle_budget","fully_feasible","feasible_shared_target_budget"],
          "length_ratio_target","Fresh saturated nuisance fits; 1000 patients/site")
    comparison=chosen[0]["length_ratio_target"]/chosen[0]["length_ratio_certified_target"]
    axis.axhline(comparison,color="#333333",linestyle="--",linewidth=1.2,label="Finite-sample target-only")
    axis.legend(loc="upper center",bbox_to_anchor=(.5,1.4),ncol=2)
    save_figure(figure,report/"fitted_cost")
    summaries=["RoCE-K：回退、校准和nuisance误差的后续实验", "",
      "完成18组高斯实验（每组2000 repeats）、18组二元分数实验（每组1000 repeats）、12组全新低维nuisance拟合（每组300 repeats）。",
      "共57600次配置内重复；nuisance预算改进为同一数据上的配对复算，不额外计作独立样本。",
      "另审计90份既有真实RoCE校准拟合，并独立做40000次target-only基准检查。所有主实验每site总样本1000。",
      "", "主要判断：",
      "1. 推荐继续保留固定错误预算下的较短保护回退。每次区间都有来源moment区间的长度上界，且不丢弃空集/回退重复。",
      "2. K依赖的target错误预算在没有共享target噪声时有益；存在共享噪声时，可能把区间变宽，不作为统一默认。",
      "3. Bernstein型target方向区域显著优于粗糙矩椭球。真实估计方差后，完整保护仍有明显效率成本。",
      "4. 实际拟合后，可计算的最坏情况偏差预算远大于精确条件偏差。共享target预算只支付一次有所改善，仍不足以获得实用借用收益。",
      "", "高斯K=2048，全部有效，无共享target底限",table(gaussian,2048,0.,"all_valid",["fixed_target","fixed_shorter","decay_target","decay_shorter"]),
      "", "高斯K=2048，强偏差，无共享target底限",table(gaussian,2048,0.,"strong",["fixed_target","fixed_shorter","decay_target","decay_shorter"]),
      "", "高斯K=2048，弱偏差，共享target SD=0.5 source SE",table(gaussian,2048,.5,"weak_boundary",["fixed_target","fixed_shorter","decay_target","decay_shorter"]),
      "", "二元分数K=2048，弱偏差，效应异质性=.3",table(binary,2048,.3,"weak_boundary",methods+["gaussian_plugin"]),
      "", "全新饱和模型拟合K=128，弱偏差，效应异质性=.3",table(fitted,128,.3,"weak_boundary",
          ["gaussian_zero_budget","gaussian_oracle_budget","bounded_oracle_budget","fully_feasible","feasible_shared_target_budget"]),
      "", "同样有限样本保护的target-only比较："]
    for row in fitted:
        if row["source_count"]==128 and row["shared_scale"]==.3 and row["pattern"]=="weak_boundary":
            summaries.append("{}: length / certified target = {:.3f}; certified target coverage = {:.3f}".format(
                row["method"],row["length_ratio_certified_target"],row["coverage_certified_target"]))
    references=json.loads((root/"target_reference_reproduction.json").read_text())
    summaries += ["", "普通target-only参考区间是渐近正态区间，不称为有限样本精确95%。独立20000 repeats/情形的结果："]
    summaries += ["heterogeneity={}: coverage={:.5f}, MCSE={:.5f}".format(r["effect_heterogeneity"],r["coverage"],r["mcse"]) for r in references]
    summaries += ["", "解释与边界：",
      "高斯部分是抽象联合残差模型；二元分数和饱和拟合部分才直接生成患者级因果数据。不能把高斯结果当作原始因果模型的完整验证。",
      "饱和模型在每个repeat中重新拟合共同OR、target PS及source校准权重；它是低维闭式诊断，不是生产版稀疏/fold-summed RoCE。",
      "source用500训练+250方差pilot+250评估；target用500训练+500评估，另构造使用全部1000患者的两折target-only参考。",
      "oracle bias budget和oracle variance仅用于定位误差来源。fully_feasible使用联合二项概率界、经验方差证书和经验Bernstein区域。",
      "主结果大多高于95%coverage，校准仍保守。较短回退改善RMSE和长度，不等于完整高维问题已经解决。",
      "90份既有RoCE拟合只有每个训练样本的一次真实留出评估；它们的覆盖计数不作为新Monte Carlo coverage报告。",
      "下一步应优先利用校准方程推导更紧的整体余项控制，并改进联合方差校准。不要把无效的零预算或未经证明的plug-in当作最终方法。",
      "文献依据：Maurer--Pontil (2009) Theorems4/10；其两侧和多方向版本使用显式联合错误概率。https://arxiv.org/abs/0907.3740"]
    (report/"REPORT_zh.txt").write_text("\n".join(summaries)+"\n")
    audit=dict(gaussian_cells=18,binary_cells=18,fitted_cells=12,reported_repeats=57600,
               archived_fits=90,reference_audit_repeats=40000,method_setting_rows=len(gaussian+binary+fitted),
               all_intervals_recomputed=True,shorter_fallback_length_cap_checked=True,paired_fitted_reanalysis_exact=True)
    (report/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    latex=[r"\begin{tabular}{llrrrr}",r"\toprule",r"Scenario & Rule & Coverage & Length/target & RMSE midpoint & RMSE source\\",r"\midrule"]
    for pattern in ["all_valid","weak_boundary","strong"]:
        for method in ["fixed_target","fixed_shorter","decay_shorter"]:
            r=next(r for r in gaussian if r["source_count"]==2048 and r["shared_scale"]==0 and r["pattern"]==pattern and r["method"]==method)
            latex.append("{} & {} & {:.3f} & {:.3f} & {:.5f} & {:.5f}\\\\".format(pattern.replace("_"," "),
                {"fixed_target":"Original","fixed_shorter":"Shorter","decay_shorter":"Decay + shorter"}[method],
                r["coverage"],r["length_ratio_target"],r["rmse_midpoint"],r["rmse_source"]))
    latex += [r"\bottomrule",r"\end{tabular}"]
    (report/"gaussian_table.tex").write_text("\n".join(latex)+"\n")
    print(json.dumps(audit))


if __name__=="__main__":main(Path(sys.argv[1]))

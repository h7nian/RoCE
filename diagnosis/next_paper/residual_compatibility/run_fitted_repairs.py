"""Fresh finite-stratum nuisance fitting and archived RoCE-fit diagnostics."""
import argparse
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import time
import numpy as np
from fitted_strata import fitted_stratum_data
from repairs import infer_repair, variance_certificate
from run_repairs import save_case


def fitted_plan():
    rows=[]
    for count in [8,32,128]:
        for pattern in ["all_valid","weak_boundary"]:
            for floor in [0.,.3]:
                rows.append(dict(cell_id=len(rows)+1,panel="fitted_strata",source_count=count,pattern=pattern,
                    shared_scale=floor,n_per_site=1000,source_training_size=500,source_pilot_size=250,
                    source_evaluation_size=250,target_training_size=500,target_evaluation_size=500,
                    valid_fraction=.75,seed_offset=6200000))
    return rows


def run_fitted(task):
    setting,repeats,output=task
    start=time.time()
    data=fitted_stratum_data(setting,repeats)
    results={}
    oracle_target=dict(data,target_sample_covariance=data["target_covariance"])
    results["gaussian_zero_budget"]=infer_repair(oracle_target,target_calibration="normal")
    results["gaussian_oracle_budget"]=infer_repair(oracle_target,target_calibration="normal",
        source_bias_budget=data["oracle_source_budget"],target_bias_budget=data["oracle_target_budget"])
    results["bounded_oracle_budget"]=infer_repair(data,target_calibration="empirical",bounded=True,
        source_bias_budget=data["oracle_source_budget"],target_bias_budget=data["oracle_target_budget"])
    results["bounded_feasible_budget"]=infer_repair(data,target_calibration="empirical",bounded=True,bias_failure=.005,
        source_bias_budget=data["feasible_source_budget"],target_bias_budget=data["feasible_target_budget"])
    lower,upper=variance_certificate(data["pilot_variance"],data["score_ranges"],250,250,.005)
    results["fully_feasible"]=infer_repair(data,target_calibration="empirical",bounded=True,
        target_alpha=.02,variance_lower=lower,variance_upper=upper,variance_failure=.005,bias_failure=.005,
        source_bias_budget=data["feasible_source_budget"],target_bias_budget=data["feasible_target_budget"])
    results["feasible_shared_target_budget"]=infer_repair(data,target_calibration="empirical",bounded=True,
        target_alpha=.02,variance_lower=lower,variance_upper=upper,variance_failure=.005,bias_failure=.005,
        source_bias_budget=data["shared_source_budget"],target_bias_budget=data["shared_target_budget"])
    metrics=save_case(output,setting,data,results,time.time()-start)
    diagnostics=dict(setting=setting,
        mean_oracle_source_budget=np.mean(data["oracle_source_budget"],axis=0).tolist(),
        mean_feasible_source_budget=np.mean(data["feasible_source_budget"],axis=0).tolist(),
        mean_oracle_target_budget=np.mean(data["oracle_target_budget"],axis=0).tolist(),
        mean_feasible_target_budget=np.mean(data["feasible_target_budget"],axis=0).tolist(),
        mean_shared_source_budget=np.mean(data["shared_source_budget"],axis=0).tolist(),
        mean_shared_target_budget=np.mean(data["shared_target_budget"],axis=0).tolist(),
        bias_budget_failure=float(np.mean(np.any(data["oracle_source_budget"]>data["feasible_source_budget"],axis=1)|
             np.any(data["oracle_target_budget"]>data["feasible_target_budget"],axis=1))),
        variance_certificate_failure=float(np.mean(np.any((data["variance"]<lower)|(data["variance"]>upper),axis=(1,2)))))
    (Path(output)/"diagnostics"/("cell_{:03d}.json".format(setting["cell_id"]))).write_text(json.dumps(diagnostics,indent=2)+"\n")
    return metrics


def archived_diagnostics(input_file,output):
    packets=json.loads(Path(input_file).read_text())
    rows=[]
    for packet in packets:
        data=dict(source=np.array(packet["source"])[None,:,:],variance=np.array(packet["variance"]),
                  target=np.array(packet["target"])[None,:],target_covariance=np.array(packet["target_covariance"]),
                  target_sample_covariance=np.array(packet["target_covariance"]),target_sample_size=packet["evaluation_sizes"][0])
        for covariance_mode in ["population","empirical"]:
            current=dict(data)
            if covariance_mode=="empirical":
                current["variance"]=np.array(packet["empirical_variance"])
                current["target_covariance"]=np.array(packet["target_empirical_covariance"])
                current["target_sample_covariance"]=current["target_covariance"]
            for budget in ["zero","oracle"]:
                result=infer_repair(current,target_calibration="normal",valid_minimum=packet["valid_minimum"],
                    source_bias_budget=packet["source_bias_budget"] if budget=="oracle" else None,
                    target_bias_budget=packet["target_bias_budget"] if budget=="oracle" else None)
                interval=result["interval"][0]
                rows.append(dict(packet["setting"],covariance_mode=covariance_mode,budget=budget,
                    lower=float(interval[0]),upper=float(interval[1]),truth=packet["truth"],
                    covers=bool(interval[0]<=packet["truth"]<=interval[1]),
                    source_budget_max=max(packet["source_bias_budget"]),
                    target_budget_max=max(packet["target_bias_budget"]),fallback=int(result["fallback"][0])))
    with (output/"archived_fit_diagnostics.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (output/"COMPLETE").write_text("90 archived fitted samples audited; not a new Monte Carlo coverage experiment.\n")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--repeats",type=int,default=300)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--archived-packets")
    parser.add_argument("--smoke",action="store_true")
    args=parser.parse_args();output=Path(args.output)
    if args.repeats<1 or not 1<=args.workers<=4:
        raise ValueError("Positive repeat count and one to four workers required")
    if not str(output).startswith("/scratch.global/zhan9381/FACE-HD/") or output.exists():
        raise ValueError("Use a new FACE-HD scratch output directory")
    output.mkdir(parents=True)
    if args.archived_packets:
        archived_diagnostics(args.archived_packets,output)
        return
    (output/"cells").mkdir();(output/"diagnostics").mkdir()
    plan=fitted_plan()
    if args.smoke:plan=plan[:1]
    hashes={file.name:hashlib.sha256(file.read_bytes()).hexdigest() for file in Path(__file__).parent.glob("*.py")}
    (output/"configuration.json").write_text(json.dumps(dict(settings=plan,repeats=args.repeats,code_sha256=hashes,
        scope="Fresh saturated cell fits, not production sparse fold-summed RoCE; feasible vs oracle budgets distinguished"),indent=2)+"\n")
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for metrics in pool.imap_unordered(run_fitted,[(row,args.repeats,str(output)) for row in plan]):
            records.extend(metrics);print("Completed fitted cell",metrics[0]["cell_id"],flush=True)
    with (output/"metrics.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=sorted(set(key for row in records for key in row)))
        writer.writeheader();writer.writerows(records)
    (output/"COMPLETE").write_text("All nuisance fits, held-out observations and repeats retained.\n")


if __name__=="__main__":main()

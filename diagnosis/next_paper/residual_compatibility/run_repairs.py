"""Fresh-seed paired tests of fallback and confidence calibration repairs."""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
import numpy as np
from scipy.stats import norm
from experiments import gaussian_data, binary_data
from repairs import infer_repair, variance_certificate


def make_plan(profile):
    rows = []
    for count in ([128, 512, 2048] if profile == "gaussian" else [32, 512, 2048]):
        for pattern in ["all_valid", "weak_boundary", "strong"]:
            for floor in ([0., .5] if profile == "gaussian" else [0., .3]):
                rows.append(dict(cell_id=len(rows)+1, panel=profile, noise=profile if profile=="gaussian" else "binary",
                    source_count=count, pattern=pattern, shared_scale=floor, n_per_site=1000,
                    source_evaluation_size=1000 if profile=="gaussian" else 750,
                    target_evaluation_size=1000, pilot_size=0 if profile=="gaussian" else 250,
                    valid_fraction=.75, heterogeneous=False, frequency_scale=2., seed_offset=5200000))
    return rows


def save_case(output, setting, data, results, elapsed):
    names = list(results)
    truth = np.broadcast_to(data["truth"], (len(data["target"]),))
    reference = data["target"] @ np.array([1.,1.,-1.])
    covariance = np.broadcast_to(data["target_covariance"], (len(truth),3,3))
    variance = np.einsum("i,bij,j->b", [1.,1.,-1.], covariance, [1.,1.,-1.])
    reference_width = 2 * norm.isf(.025) * np.sqrt(variance)
    reference = data.get("reference_point",reference)
    reference_width = data.get("reference_width",reference_width)
    records = []
    for name in names:
        result = results[name]
        lower, upper = result["interval"].T
        coverage = (lower <= truth) & (truth <= upper)
        length = upper-lower
        errors = result["point_midpoint"] - truth
        squared = errors**2
        fallback = result["fallback"] != 0
        records.append(dict(setting, method=name, repeats=len(truth), coverage=float(coverage.mean()),
            coverage_mcse=float(np.sqrt(coverage.mean()*(1-coverage.mean())/len(truth))),
            mean_length=float(length.mean()), length_ratio_target=float(length.mean()/reference_width.mean()),
            length_ratio_certified_target=float(length.mean()/np.mean(data["certified_reference_width"])) if "certified_reference_width" in data else None,
            rmse_midpoint=float(np.sqrt(squared.mean())), rmse_source=float(np.sqrt(np.mean((result["point_source"]-truth)**2))),
            rmse_target=float(np.sqrt(np.mean((reference-truth)**2))), bias_midpoint=float(errors.mean()),
            coverage_target=float(np.mean(abs(reference-truth)<=reference_width/2)),
            coverage_certified_target=float(np.mean((data["certified_reference_interval"][:,0]<=truth)&
                (truth<=data["certified_reference_interval"][:,1]))) if "certified_reference_interval" in data else None,
            fallback_rate=float(fallback.mean()), target_fallback_rate=float(np.mean(result["fallback"]==1)),
            source_fallback_rate=float(np.mean(result["fallback"]==2)),
            fallback_mse_share=float(squared[fallback].sum()/max(1e-30,squared.sum())),
            mean_source_bias_budget=float(np.mean(result["source_bias_budget"])),
            mean_target_bias_budget=float(np.mean(result["target_bias_budget"][:,1:])),
            alpha_target=result["alpha_target"], alpha_source_arm=result["alpha_source_arm"]))
    path = Path(output)/"cells"/("cell_{:03d}".format(setting["cell_id"]))
    np.savez_compressed(str(path)+".npz", methods=np.array(names), truth=truth,
        interval=np.stack([results[name]["interval"] for name in names],axis=1),
        point_midpoint=np.column_stack([results[name]["point_midpoint"] for name in names]),
        point_source=np.column_stack([results[name]["point_source"] for name in names]),
        fallback=np.column_stack([results[name]["fallback"] for name in names]),
        length_envelope=np.column_stack([results[name]["length_envelope"] for name in names]),
        target_reference_width=reference_width,target_reference_point=reference,
        certified_reference_width=data.get("certified_reference_width",np.full(len(truth),np.nan)))
    path.with_suffix(".json").write_text(json.dumps(dict(setting=setting,metrics=records,runtime_seconds=elapsed),indent=2)+"\n")
    return records


def run_cell(task):
    setting, repeats, output = task
    start = time.time()
    data = gaussian_data(setting,repeats) if setting["noise"]=="gaussian" else binary_data(setting,repeats)
    results = {}
    if setting["noise"] == "gaussian":
        for allocation in ["fixed", "decay"]:
            alpha = .025 if allocation=="fixed" else .025/math.sqrt(setting["source_count"])
            for fallback in ["target","shorter"]:
                results[allocation+"_"+fallback] = infer_repair(data,target_alpha=alpha,fallback_rule=fallback)
    else:
        for calibration in ["moment","bernstein","empirical"]:
            results["known_variance_"+calibration] = infer_repair(data,bounded=True,target_calibration=calibration)
        results["gaussian_plugin"] = infer_repair(data,target_calibration="normal",
            variance_lower=data["pilot_variance"])
        lower, upper = variance_certificate(data["pilot_variance"],data["score_ranges"],
            setting["pilot_size"],setting["source_evaluation_size"],.005)
        results["estimated_variance_empirical"] = infer_repair(data,bounded=True,target_calibration="empirical",
            variance_lower=lower,variance_upper=upper,variance_failure=.005)
    return save_case(output,setting,data,results,time.time()-start)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--profile",choices=["gaussian","binary"],required=True)
    parser.add_argument("--repeats",type=int,default=1000)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--smoke",action="store_true")
    args=parser.parse_args()
    output=Path(args.output)
    if not str(output).startswith("/scratch.global/zhan9381/FACE-HD/") or output.exists():
        raise ValueError("Use a new FACE-HD scratch directory")
    if args.repeats < 1 or not 1<=args.workers<=4:
        raise ValueError("Invalid repeat count or worker count")
    output.mkdir(parents=True);(output/"cells").mkdir()
    plan=make_plan(args.profile)
    if args.smoke:plan=plan[:2]
    hashes={file.name:hashlib.sha256(file.read_bytes()).hexdigest() for file in Path(__file__).parent.glob("*.py")}
    (output/"configuration.json").write_text(json.dumps(dict(settings=plan,repeats=args.repeats,
        code_sha256=hashes,scope="Oracle residual scores; no fitted sparse-RoCE coverage claim"),indent=2)+"\n")
    records=[]
    with multiprocessing.Pool(args.workers) as pool:
        for metrics in pool.imap_unordered(run_cell,[(row,args.repeats,str(output)) for row in plan]):
            records.extend(metrics)
            print("Completed",metrics[0]["cell_id"],args.profile,"K",metrics[0]["source_count"],flush=True)
    fields=sorted(set(key for row in records for key in row))
    with (output/"metrics.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(records)
    (output/"COMPLETE").write_text("All prespecified repeats retained.\n")


if __name__=="__main__":main()

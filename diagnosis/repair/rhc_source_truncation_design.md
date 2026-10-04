# RHC source truncation: controlled finite-sample diagnosis

The user authorized truncation experiments to investigate the original min
Private-target result. Keep all5,039 patients,61 shared features,Private target,
three sources,Two-layer/one-round,joint TATE,100lambda grids,min nuisance CV and
aggregation cutoff2. The calibrated target fitting/evaluation radius stays5.

Prespecify source fitting AND evaluation radii2,3,4,5, corresponding to upper
merged-weight caps exp(radius):7.39,20.09,54.60,148.41. This is symmetric log-tilt
clipping, also raising the lower cap exp(-radius); it is not upper-tail-only
trimming. Refit the entire source calibration program at each radius. Do not
change only an evaluation score or rescale its standard error. Initial target
OR messages and both target references must remain exactly fixed within each
partition. Baselines are unchanged; reuse their matched complete results.

Use the historical default partition and seeds101,202,303:16 fits. Radius5
reproduces the saved reference in each partition before the new candidates are
interpreted. This explicit reference also checks the new driver argument. Each
fit is a separate4CPU/8GB/1h Slurm job,checkpoint/requeue,preempt/public/saffo-2tb.
The checked v2 estimator library is unchanged; only driver configuration differs.

Review all prescribed radii and partitions: point/SE/CI,source-arm ESS,largest
record variance share,clipped fraction,ordinary balance,source/target candidate
gaps and numerical convergence. Lower SE or significance is not evidence of
smaller bias. Strong clipping can introduce bias and worsen moment balance;
real data cannot identify true bias/coverage. The theorem's population margin
conditions do not follow from choosing a smaller radius.

No patients are discarded and no primary/default/manuscript result is replaced
by this experiment. Preserve radius5 and all prior studies. Store results under
the FACE-HD scratch real_data/rhc hierarchy. Any further confirmation must use
fresh prespecified partitions and retain every prescribed result.

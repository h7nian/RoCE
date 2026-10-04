# RHC truncation: known-truth model-based sensitivity

This is an auxiliary experiment, not a replacement of the observed RHC analysis
or an estimate of its actual causal bias. Freeze the population generator before
fitting any replicate. Retain Private as target, the four observed site sizes,
61 shared features, original preprocessing, min nuisance CV with 100 lambdas,
two-layer/one-round joint TATE and aggregation cutoff 2. Compare source fitting
and evaluation log radii 2, 3 and 5; target radius stays 5. Keep every replicate.

Use a finite support of the distinct observed feature profiles, retaining their
empirical frequencies. The hypothetical target feature law is 95% empirical
Private plus 5% pooled empirical. This explicit smoothing supplies common
support and changes the hypothetical target; it does not change the real cohort.
Fit the population templates once, using a recorded seed and the original min
rule. Initial outcome templates use full-target arm-specific logistic lasso.
Propensity templates use each full site's logistic lasso. Preserve their mean
probability while scaling slopes and adjusting the intercept so that absolute
logits do not exceed 2. No simulated outcome is used to choose the template.

Four prespecified scenarios:

1. **O_case_mix:** common logistic outcome template; source feature laws are
   90% their empirical site distribution plus 10% pooled. The outcome working
   model is correct and the source joint-weight model is generally misspecified.
2. **W_overlap:** source joint laws are exponential tilts of the hypothetical
   target law, with maximum absolute true log merged weight 1.4. Both source
   arms have a correct linear exponential weight model. The outcome template
   has an omitted bounded age-by-APS interaction.
3. **W_tail:** the same construction allows maximum absolute true log weight
   2.4. Radius 2 can truncate the true weights; radii 3 and 5 contain them.
4. **W_departure:** W_tail, with a +1 treated outcome log-odds shift in Medicaid.
   Its control outcome is unchanged; the target truth is unchanged. This is
   a source transportability stress test, not an all-sources-valid setting.

For W scenarios, preserve each outcome arm's target mean while bounding the
linear logit template by 1.4, then add opposite +/-0.75 times the centered
interaction tanh(age) tanh(APS). Source tilt directions come from their full-data
initial merged-weight fits at radius 5. Scale each site's two directions together
to satisfy the stated true-weight bound while preserving the observed treatment
fraction. Explicitly, P_s(X=x,A=a) is proportional to f_t(x) exp(s phi(x)'gamma_a),
normalized within each arm and multiplied by its arm probability. Consequently
q_sa(x)=f_t(x)/P_s(X=x,A=a) is exactly exponential in the same features.

Each replicate samples independent X, A and Y from these fixed laws. Repeated
profiles are independent new patients, with new treatments and outcomes, not
copies of one observed outcome. The population TATE is the exact finite sum of
m1-m0 under f_t; it is fixed across replicates. An oracle target AIPW comparison
has an independently computed exact sampling variance as a generator check.

Start with 20 replicates per scenario, one job per scenario/replicate and all
three radii paired within that job. First release only one replicate per scenario;
require finite results, matching target fits across radii and successful resume.
Then release the remaining prescribed replicates. Report bias, RMSE, empirical
SD, average SE, coverage with Monte Carlo uncertainty, and paired squared-error
differences. Twenty replicates are a pilot, not adequate evidence of 95% coverage;
expand a frozen, validated design before making precise coverage claims.

The first frozen pilot exposed an initial-weight CV grid failure in W_overlap.
Retain that pilot and its failed replicate. A revised implementation may retry
an entirely failed initial grid once, on the identical CV folds, appending
larger penalties up to a bound based on feature extrema. It retains the original
100 candidates, min selection and convergence criteria; successful old paths
are unchanged. Freeze this implementation separately and record every retry.
The population template, replication seeds and truncation choices remain fixed.

The finite-support generator has bounded covariates and exact specified nuisance
truths; it does not establish all nuisance-rate or uniform inference assumptions.
Fitted-model simulations cannot assess unmeasured confounding or prove actual
RHC causal bias is small. Keep adverse results. Never select a cutoff because
its real-data confidence interval excludes zero, and do not remove influential
patients from the primary cohort.

Rationale: weight truncation trades variance against residual bias (Cole and
Hernan, 2008, <https://pmc.ncbi.nlm.nih.gov/articles/PMC2732954/>). Model-based
positivity diagnostics can be optimistic because they assume fitted models
(Petersen et al., 2012, <https://pmc.ncbi.nlm.nih.gov/articles/PMC4107929/>).
These experiments diagnose a proposed numerical/statistical modification; they
do not bias-correct the observed RHC effect.

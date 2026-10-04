# RHC covariate and weight diagnosis

The user authorized deeper diagnosis and changes to confounder specifications.
The two new profiles below are exploratory responses to the initial RHC
results. They do not replace the historical specification or establish a
confirmatory clinical finding. All profiles and unfavorable findings remain.

## Findings from saved fits

The primary risk difference is.07895 with SE.04511. Saved two-arm influence
functions and exact deployed source weights reproduce the reported variance
and source corrections. Medicaid contributes80.0% of total variance. A single
Medicaid treated record contributes63.0% of total variance and has weight148.4,
the exp(5) cap. The Medicaid treated arm has193 records but a descriptive Kish
weight ESS of25.0. The other source/arm ESS values range from37.8 to91.6 despite
490-947 raw records. This indicates strongly concentrated weights, not evidence
that final optimizers failed: final-fit nonconvergence counts are zero.

The main Medicaid record has urine missing, white-cell count0 and primary
MOSF with malignancy. Zero white cells may represent an extreme clinical value;
neither that record nor its outcome is deleted or recoded as an error. The
diagnosis concerns extrapolation, residual weighting and sampling instability.
The two candidate mean/covariance identities passed, so the observed larger
SE is not explained by having simply omitted the arm covariance.

The primary clipped fraction is only6/3341 source records across arms. Small
clipping fractions can still conceal large influence concentration. Ordinary
weighted SMDs reach about.6. These descriptive SMDs are not the derivative-
weighted moment constraints of our calibration loss; a large SMD does not by
itself demonstrate an implementation error in that loss.

## Confounder timing and literature

[Connors et al. (1996), pp.890-891](https://ekroc.weebly.com/uploads/2/1/6/3/21633182/connorsetal1996.pdf)
describe RHC during the initial24h and physiological summaries over that day.
This does not establish every physiological measurement preceded treatment.
The original propensity-model description explicitly identifies ADL/DASI as
function two weeks before admission, despite later interviews. DASI is
therefore retained; collection date alone is not a reason to label it a mediator.
Important severity, demographic, disease and prior-history covariates remain.

The [public variable catalogue](https://hbiostat.org/data/repo/crhc) documents
urine-output missingness. In the retained cohort2651/5039 urine values are
missing. A missingness indicator is a sensitivity analysis, not a proof that
missing-data bias or temporal ambiguity is solved.

[Hirano and Imbens (2001)](https://link.springer.com/article/10.1023/A:1020371312283)
is in Health Services and Outcomes Research Methodology, not Biostatistics.
Our historical61-feature insurance-stratified matrix is not claimed to be an
exact reconstruction of their72-term published design.

## Frozen covariate profiles before the new fits

All5039 patients, treatment/outcome definitions and the Private target remain.
Age, DASI, prognosis, APACHE/coma, physiological and comorbidity information
remain in the models. OR and joint weights always share the same feature map.

1. `log_missing` (62 features): log1p urine output, creatinine, bilirubin and
   white-cell count after the existing median imputation; add `urin1_missing`.
   Other covariates are unchanged. The transformations preserve the observed
   numerical information and modify the working functional forms.
2. `grouped_log_missing` (60 features): additionally combine primary Colon
   Cancer/Lung Cancer into Solid Cancer and combine the trauma/orthopedic flags
   into `injury_diagnosis`. Cancer status and prior malignancy remain separate.
   Grouping coarsens information and may leave residual confounding; it is a
   sensitivity profile, not an established sufficient adjustment set. It removes
   empty source-arm categories in the full-cohort marginal table, but sparse
   fold-specific cells and poor joint support can remain.

No target trimming, individual deletion or outcome-dependent variable selection
is applied. Fully observed age is not removed to conceal differences between
insurance populations. Pooled covariate preprocessing remains a documented
limitation of this retrospective federated emulation.

## Evaluation and execution

Run the same12 tasks for each new profile: primary, four paired baselines,
source-standard/ordinary-target calibrations, CV1se, radius12 and three fold
seeds. The primary uses Two-layer, one-round, joint TATE, cutoff2,100 lambdas,
PN tolerance1e-10, and fitting/inference radius5. Reuse fitted models for the
three aggregation modes at cutoffs1/2/3. Compare point estimates, CIs, two-arm
weights, support flags, source/arm ESS, site variance contributions and maximum
individual variance share. Do not rank profiles by statistical significance or
claim smaller SE establishes smaller bias.

The original jobs took5-9min with about0.8-1.2GB reported MaxRSS. The new jobs
request4 CPUs,8GB and1h, with checkpoints/requeue and the requested public plus
saffo-2tb routes. Exclude acl45/acl47 as well as the previously reviewed faulty
nodes. New results live under `real_data/rhc/covariate_sensitivity_v1/`.

The historical ordinary-target ablation's driver-control error passed38 checks
and was separately resubmitted as job2201380 under `ordinary_target_fix_v1/`.
The old failed attempt remains. Subsequent covariate interface checks passed59
assertions, including exact historical cohort equality, identical patient/
treatment/outcome rows, outcome-invariant coding and shared OR/weight features.

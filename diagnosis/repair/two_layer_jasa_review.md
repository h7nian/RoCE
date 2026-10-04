# Two-layer review of the current JASA manuscript

Read-only snapshot: Overleaf project67abae7087d0ce709abdbb70, revision
`d57e90e5ba616a404e35d04da812fb63ef244b0d`, file `JASA/JASA.tex`.
The snapshot is under scratch `implementation/r9/jasa_review_20260921/`.
No manuscript changes were pushed or written into the user's checkout.

## Verdict and correction of the earlier interpretation

The JASA construction is reasonable as a TWO-LAYER setting: k1 is outer
evaluation and k2 is calibration. Initial fitting on the complement of those
two folds does not introduce a third fold index or a third validation loop.

The manuscript explicitly retains fold-specific eta^(-k1), learned only in
T_(-k1)=D excluding outer fold k1. This is stated in active text at lines662
and677, and again in the expanded algorithm at1156. It is NOT the global-eta
pooled-outer-score proposal mentioned earlier in the conversation. The layer
ablation should follow the JASA fold-specific definition unless the user
subsequently asks for a separate global-eta comparison.

## What is already specified

- Initial models for calibration block k2 exclude both k1 and k2.
- The final outer nuisance models minimize the sum/average of all k2 losses;
  parameters from separate fits are not averaged (lines1033-1063).
- Eta and any aggregation-cutoff selection exclude k1 (lines662,677).
- Final estimates use the untouched outer fold, then combine outer results.
- K and K_f are fixed in the current theory, while nuisance dimension p may
  increase (lines698-716). Small site count is not the same as small p.

## The missing implementation detail

Line1156 says the training-sample moments use the same secondary-fold
cross-fitting scheme, but it does not unambiguously name the nuisance model
used for each eta-learning score. The text must distinguish:

1. Initial plug-in scores, whose fitted models exclude k1 and k2.
2. Scores from the final outer-calibrated model, evaluated inside T_(-k1).
3. Independently validated final-calibrated scores, whose training would
   require an additional calibration index k3.

The first changes the candidate program being judged. The third is the
current three-layer construction. For the advisor's two-layer ablation, the
recommended complete definition is the SECOND option: evaluate the final
outer-calibrated source and target models inside T_(-k1) for weight learning,
with no additional independent eta-validation fold. These are training-sample
scores; do not label them independent nuisance-validation scores.

The source evaluation summaries immediately preceding line1156 are defined
on outer fold k1. Weight learning instead needs separate source TRAINING
correction means, variances and cross-arm moments. Their data labels must
exclude k1. They can be sent in the existing batch with separately labeled
outer-evaluation summaries; eta fitting must not consume the latter.

The target summaries also have two different roles: the initial weight fit
uses the target feature mean excluding k1,k2; the calibrated weight loss uses
the derivative moment on target calibration fold k2, evaluated with its
out-of-k1,k2 initial OR. The expanded description makes this distinction;
the main paragraph's shorter wording should not obscure it.

## A coherent fixed-K proof route

For valid sources, write D_(j,k1) as the source/target difference evaluated
on the outer fold. Under the calibrated score expansion, D_(j,k1)=O_p(n^-1/2).
If eta_hat^(-k1)->eta* and incompatible sources are excluded with probability
tending to one, then for fixed K and K_f,

    sum_j (eta_hat_j^(-k1)-eta*_j) D_(j,k1) = o_p(n^-1/2).

Only consistency of eta is needed here, not a root-n rate for eta. The final
outer data are independent of both fitted nuisances and eta because both are
functions of T_(-k1). Eta need not have another independent validation split
from the nuisance fits inside that training sample.

However, consistency of the eta-training moments must be established for the
actual chosen scores. If final nuisance fits and eta moments reuse T_(-k1),
one cannot simply invoke an independent-evaluation-fold law of large numbers.
A sufficient route is empirical prediction consistency for the training
scores, plus uniform empirical-gradient control for their means. For example,
with population score orthogonality, an empirical gradient sup-norm of order
sqrt(log(p)/n), nuisance L1 error of order s*sqrt(log(p)/n), and a controlled
quadratic remainder, the extra mean error is O_p(s*log(p)/n). The usual
s*log(p)=o(sqrt(n)) condition then makes it negligible. Covariance consistency
can follow from empirical L2 score consistency and Cauchy-Schwarz. These are
conditions to prove, not an automatic consequence of counting two layers.

The current manuscript's oracle proof additionally uses a growing cutoff
1/lambda_N, source separation and a stable oracle minimizer. A fixed cutoff2
does not itself prove eta consistency. Keep finite-sample variance checks
separate from that asymptotic statement.

## Controlled two-versus-three-layer experiment

The implemented public argument is `crossfit_layers=2L` or `3L`. The checked
library is frozen under scratch `implementation/r9/layer_checks_v1/`; the
paired experiments are under `r9/layer_validation_v2/`. All p10/p100/p200
campaigns have completed450 results per layer and dimension,2700 total.
The final report is `r9/layer_review_p200_complete_v1/`; existing simulation
campaigns retain their frozen implementations. See
`layer_validation.md` for the design and timestamped scratch reports for results.

Both variants retain the SAME outer fold partition, final outer calibration
losses, nuisance lambda grid, solver tolerance, aggregation objective, cutoff,
and target anchor. Preserve the existing three-layer default until evidence
supports a documented change. One/two communication rounds remain a separate
argument.

| Item, 1000/site and10 folds | Two-layer proposal | Existing three-layer |
|---|---|---|
| Final outer nuisance training | 900 rows; each initial block fit800 | Same |
| Models used to learn eta | Reuse those final outer models | Independently trained800-row models, each initial fit700 |
| Observations supplying eta summaries | The900 outer-training rows | The same900 rows, through rotating inner validation |
| Final evaluation | The100 outer-held-out rows | Same |
| Eta | One fold-specific vector per arm/mode | Same fold-specific form |

Use the same moment aggregation/centering convention for both variants.
Require final outer nuisance coefficients, target-only estimates and raw
source-candidate estimates/variances to match. Changes in aggregated means
and TATE should arise from eta learning, not changed nuisance definitions.

Report eta distributions, aggregated mu0/mu1/TATE, their covariance,
estimated fixed-weight IF variance, eta-sensitivity-adjusted variance,
empirical variance across repeats, SE/SD ratios, coverage, bias, RMSE and time.
The existing eta sensitivity holds nuisance fits and active sets fixed; it
must not be described as a full-refit bootstrap or universally valid under
selection kinks. Empirical repeated-simulation variance is the primary check.

Vary source count K separately from covariate dimension p. Fixed K=2,4,6 is
the immediate comparison; growing-K settings belong to the subsequent study.
Three layers may help control a high-dimensional weight learner, but do not
by themselves prove growing-K validity or fix weak-source selection.

## Suggested clarification for the expanded algorithm

The following is a proposed completion of the eta step, not a quotation from
the manuscript or a claim that its existing text already specifies this choice:

> For each outer fold k1, we first obtain the final calibrated nuisance fits
> using the secondary-fold loss aggregation described above. We reuse these
> final fits to evaluate estimating-function contributions on the outer
> training sample T_(-k1). Each source supplies separately labeled training
> correction means, variances, and cross-arm moments, in addition to its
> outer-evaluation summaries. The target combines the training summaries with
> its own training-sample contributions and solves for eta^(-k1). No observation
> in outer fold k1 enters any statistic used to learn these weights. The
> secondary folds organize the calibration losses and training statistics;
> they are not independent validation folds for the final calibrated fits.

For fixed K, an accompanying result must establish consistency of these
training-sample means and covariance matrices under the fitted nuisance
procedure. The final outer holdout alone does not prove that intermediate
statement. Optional tuning of the aggregation cutoff would need its own
explicit definition; the controlled comparison fixes cutoff2.

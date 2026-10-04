# A common-outcome candidate to test

This is a research construction, not a change to the current paper or a novelty
claim. It may avoid the extra training-mean replacement term in the translated
source candidate. Evaluate it on the same saved fitted data before choosing a
research default.

For arm a, let m_beta(x)=expit(phi(x)' beta). Fit a common beta at the target
using inverse target propensity weights I(A=a)/pi_a(X). Each source retains
its own calibrated merged weight q_j. Consider

    S_j(beta,gamma_j) = E_t m_beta
                      + E_j I(A=a) q_gamma_j (Y-m_beta).

Its sample version uses the same fitted m_beta in the target prediction and
every source residual. Thus all source candidates share one target evaluation
component without adding a training-sample shift. The target anchor can remain
the separately calibrated Hou estimator.

## Population orthogonality

Assume a valid source, shared OR/weight features phi, and the usual
transportability/identification conditions. For the exponential weight model,

    partial_beta S_j = E_t m'_beta phi - E_j I_a q_j m'_beta phi,
    partial_gamma S_j = -E_j I_a q_j (Y-m_beta) phi.

These derivatives vanish under either of the following branches.

**Correct joint weights and correct target propensity.** The first derivative
vanishes for every beta. Correct target IPW makes beta* the unpenalized target
population logistic projection, so E_t (m_true-m_beta*) phi=0. Transporting the
second derivative with the correct q_j then also gives zero. The common OR
need not be correctly specified. Ordinary target arm regression or Hou's
odds-weighted OR generally does not supply this target projection.

**Correct common OR.** Any bounded positive target propensity working weight
has beta_true as the population OR minimizer. The second derivative is zero
by the conditional outcome mean. If the source calibration balances
m'_beta_true phi, the first derivative is also zero even if q_j is misspecified.
Existing calibration uses an initial OR derivative; its probability limit must
therefore equal beta_true in this branch. Ordinary target initial OR has that
property for a valid correctly specified common outcome model.

When q is clipped, the gamma derivative also has the clipping derivative.
The correct-OR branch still makes it zero. The correct-weight branch requires
the claimed true exponential weight to lie within the clipping range; the
population calibration equations and required curvature must be checked.

This is a population lemma. It does not prove that CV-selected penalties,
estimated propensity weights, empirical calibration or growing-K nuisance
remainders satisfy the rates required for confidence intervals.

In particular, correct source joint weights alone do not establish the second
derivative condition when both the target propensity and common OR are
misspecified. The target projection requirement is substantive. Do not assume
this construction inherits every robustness configuration of the original
source-specific-OR method.

## Difference from the translated-source construction

For a source-specific fitted m_j and the same fitted q_j,

    mu_common - mu_source_specific
      = (P_target_eval - P_source_eval I_a q_j)(m_common-m_j).

There is no P_target_train prediction shift. Population orthogonality supplies
a route to a second-order difference for valid sources when the relevant
nuisance limits agree. In finite samples this candidate can still have
substantial conditional nuisance bias; neither shared target predictions nor
independent evaluation removes that bias automatically.

For invalid sources, the residual mean can be nonzero and the two-arm private
covariance must be retained. A source-count confidence procedure must cover
those arbitrary invalid means and account for a uniform bound/rate on the
remaining valid-candidate biases. Independent source residuals conditional on
training and a shared target prediction are structural conveniences, not a
complete inference theorem.

## Computation and communication

The target can send its final common OR coefficients together with the
existing calibration messages. Each source evaluates its own q_j(Y-m_common)
and returns residual moments. This does not require another communication
exchange. A dedicated implementation might omit final source OR fitting, but
the initial ORs needed by the weight calibration remain. That simplification
must be validated separately; the first comparison reuses existing q fits and
changes only the outcome model in the evaluation score.

The first paired conditional-population audit is under scratch
implementation/next_paper/v9/honest_population_audit_v3/. On its one saved C3
training sample, the source-s2 TATE conditional bias changes from -.0157 for
translated_source to .0030 for common_outcome. On the C2 sample, it instead
increases from .0036 to .0080. These are single-training-sample diagnostics,
not a coverage or superiority result. Both constructions and the original
source-specific reference remain recorded. No dedicated common-outcome fitting
program or calibrated confidence interval has been released.

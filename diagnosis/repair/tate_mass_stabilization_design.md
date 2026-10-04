# Prespecified exploratory aggregation stabilization

Motivation: tighter source clipping can prevent population calibration moments
from being attainable. In the frozen O_case_mix law, a Medicaid treated income
moment already exceeds what any weight bounded by exp(3) can supply. Keep that
adverse result. Do not equate smaller reported SE with justified inference.

Evaluate one alternative without changing the nuisance losses, Private target,
min rule, two-layer folds, joint-TATE optimizer or existing radius5 fits. This
is a new optional estimator variant, not a claim that the original method has
already passed additional validation. No manuscript/default changes are made.

For each outer fold, a source returns its mean held-out merged weight including
the treatment indicator, d_ja=mean_s[I(A=a) q_ja(X)]. These are computable from
evaluation X,A and the fitted weights, without evaluation outcomes. Define

\[
g_k=1/\max(1,\max_{j,a}d_{ja,k}).
\]

Multiply both arm-specific source weight vectors by the same g_k. Scale each
source *increment relative to the target anchor*, including its target prediction
term, rather than modifying only its residual correction. The threshold1 is
fixed before examining the new variant. Do not tune it for significance. Extra
communication consists of these scalar sums; no additional round is required.

With unequal fold sizes, use the existing sitewise pseudovalue assembly. A
naive average of scalar fold estimates is not an exact substitute. Scale each
outer fold's projection of the existing eta-learning derivative by g_k as well.

First-order justification requires the original valid-candidate/inference
conditions. If the nuisance limits are shared across folds, each mass converges
to a finite constant, so g_k converges to one common g_* in (0,1]. For valid
selected sources their mean differences from the anchor are O_p(n^-1/2).
Replacing g_k by g_* consequently costs o_p(n^-1/2); this argument also requires
invalid-source contributions to be negligible under the original selection
conditions. With correct merged weights all d_ja limits equal1 and g_*=1, so
the modification is first-order equivalent to the original estimator.

If X is the anchor's limiting error and D the original aggregate-minus-anchor
error, the stabilized error is X+g_*D. Variance is convex in g. Hence an original
variance no larger than the anchor's remains no larger after this common
contraction. This does not assert that contraction always improves the original
variance or that estimated finite-sample weights have reached their limits.

The research helper uses the corresponding first-order score, including the
rescaled existing eta sensitivity. It omits the mass-factor derivative because
its population multiplier is zero under those conditions. That omission is not
a finite-sample coverage guarantee, nor a remedy for nonnegligible invalid
sources, nuisance remainders or selection kinks. Validate empirically before
calling its standard error reliable or presenting it as a main result.

Required checks: scale1 exactly recovers saved points and both SEs; scale0
recovers the calibrated target anchor; constant scales reproduce a linear
combination of the saved influence vectors; deployed weight sums reconstruct
the original source residuals. Then inspect original RHC and use the already
prescribed known-truth repeats, paired with unchanged radius5 fits. Preserve
all truncation results and any unfavorable stabilization results.

The original-partition diagnostic passed the identity, anchor, residual and
influence-linearity checks. The mass-contracted estimate is4.022pp with
first-order SE2.527pp, versus target-only SE2.456pp. This did not supply the
desired precision improvement and would add an inference condition to the
application. It was therefore not adopted or advanced as a main result. Its
Monte Carlo coverage has not been validated. The simpler radius3 procedure
has its own complete prescribed200-repeat validation; all mass-contraction
code and unfavorable diagnostic results remain as research artifacts.

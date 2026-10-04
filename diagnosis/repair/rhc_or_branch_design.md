# Additional OR-correct branch check, without replacing the stress test

The original O_case_mix law remains a severe overlap stress test. Its known
radius3 moment-feasibility violation means it cannot alone verify the intended
OR-correct inference branch under all calibration assumptions.

Add one prespecified law, O_moderate_mix. Retain the target law, original full-
target logistic outcome coefficients, bounded site propensities, site sizes,
feature map, all fitting controls and radii2/3/5. Replace each source feature
law by 90% of the target law plus10% of its original empirical site law. No new
outcomes or previous Monte Carlo estimates determine these probabilities.

The target law already puts at least0.05/5039 on every profile; the smallest
source has647 observations. For source mixture fraction0.10 and propensities
bounded by logistic(+/-2), true merged weights lie between approximately0.069
and9.32, strictly inside [exp(-3),exp(3)]. The exact finite-population check
verifies these bounds separately for both arms and every source, and verifies
the true q balances initial and outcome-derivative feature moments. Full design
rank and an interior feasible weighting law establish the relevant convex
moment problem's feasibility. A nonzero projection residual of log(q) onto the
shared linear design verifies that the working joint-weight model is wrong.
The outcome model is exactly the shared logistic working model.

Run200 repeats with the original seed formulas and IDs1–200, paired across
radii. These are new source samples under a new law; the target samples match
the existing OR-correct law, so comparisons are paired. Preserve every severe-
shift result. Do not pool both laws into one coverage statistic. This check
supports the clean OR branch; it does not prove real RHC bias or eliminate
finite-sample rate conditions.

# Independent covariance pilot for the fitted-score diagnostic

This is a next-paper diagnostic, separate from the current paper's two- and
three-layer cross-fitting comparison. It does not add a layer to either
production algorithm or change its default.

## Sample roles and implementation

Reuse the previously fitted honest models: each site has 1,000 observations,
500 of which trained the models. Split the original 500 held-out observations
into 250 covariance-pilot rows and 250 final mean-evaluation rows. The split
uses only the original observation indices. It does not balance treatment,
outcome or fitted predictions. Under the diagnostic's iid sampling assumption,
the two held-out groups are independent conditional on training.

`evaluate_honest_candidates(..., evaluation_indices=NULL)` preserves the
original full-held-out evaluation exactly. An explicit list must contain each
site once, with at least two distinct held-out indices per site. Training rows
are rejected, and the existing training-data hashes are still checked.

`split_honest_evaluation()` supplies the two disjoint index lists.
`rescale_honest_covariance()` converts pilot score covariances to the sample
sizes used for final means. For site j, the existing moment convention is

    covariance_pilot_j = centered_crossproduct_j / n_cov_j^2,
    covariance_for_final_j = covariance_pilot_j * n_cov_j / n_eval_j.

The common target block and each private source block must be scaled separately.
A single global multiplier is wrong when site sample-size ratios differ.
Cross-arm covariances are retained. No ridge, clipping, or Bessel correction
is introduced. In particular, this retains the existing finite-sample downward
factor (n_cov_j-1)/n_cov_j in the empirical covariance expectation.

## Checks and completed diagnostic

The focused check script uses a saved C2, K=8 fit. Ninety assertions passed:
exact default-packet compatibility, disjointness, index validation, unchanged
random state, mutation-based independence checks, unchanged training shifts,
training-data rejection and direct covariance reconstruction with unequal site
sizes. Artifacts are in scratch
`implementation/next_paper/v11/honest_covariance_pilot_checks_v1/`.

`review_honest_covariance_pilot.R` then reused all 90 prespecified fitted cases
from the earlier p=4, C1–C3, K=2/8/16, rho=0/.5 diagnostic. Both candidate
constructions exactly reproduced their original packets: 180 parity checks.
Original population covariance components were rescaled to the new evaluation
sizes; no additional nuisance fitting or quadrature approximation was needed.
The full report is in scratch
`implementation/next_paper/v11/honest_covariance_pilot_review_v1/`.

Independent-pilot oracle-valid GLS had reported/actual conditional variance
ratios .865–1.057 (median .978). The realized relative covariance error
`max |eigen(Sigma_hat^(-1/2) Sigma Sigma_hat^(-1/2))-1|` ranged from .078 to
.584. These errors matter even though the scalar GLS variance ratios look
closer to one: the dispersion statistic uses more covariance directions.

These are five training draws per cell, with oracle validity labels and
population moments. They are not a coverage experiment or a high-probability
covariance bound. For pilot-based weights, conditioning on training and pilot
gives the stated moments of the final estimate. The same-evaluation comparison
only reports moments with its realized weights held fixed; its weights depend
on the final observations.

## Remaining inference work

The independent split supplies the independence premise of the
[block-score approximation route](honest_block_score_clt.md). It does not
supply the covariance concentration event or valid-score bias envelopes needed
by that argument. The empirical covariance errors should be used to assess
such bounds before claiming feasible intervals. The smaller evaluation sample
also increases noise: smaller standardized drift after splitting is not
evidence that nuisance bias was reduced.

Reusing training rows for covariance could avoid this sample cost, but would
require a different proof. Do not silently substitute that design for the
checked independent-pilot construction.

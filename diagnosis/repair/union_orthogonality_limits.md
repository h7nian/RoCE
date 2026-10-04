# Initial/final nuisance limits under the union model

This concerns the uniform initial/final-limit alignment statement in the
reviewed JASA revision d57e90e5, not a claim about subsequent manuscript edits.
The two-layer data-splitting definition remains reasonable.

For a transportable source and one arm, write the population score as

    Psi(alpha,gamma) = E_t m_alpha(X)
                       + E_s[I(A=a) q_gamma(X){Y-m_alpha(X)}].

With a common working feature map phi, its derivatives involve

    d_alpha Psi = E_t[h_alpha phi] - E_s[I q_gamma h_alpha phi],
    d_gamma Psi = E_s[I (d_gamma q_gamma){Y-m_alpha}].

In the current population check, truncation is inactive, and q is an
exponential tilting model. The sign convention for log q does not affect the
zero-gradient statements.

## Correct joint-weight model

Both the initial and calibrated weight limits can equal q_true. This is
because the true density/treatment ratio balances target moments for every
integrable multiplier, including the initial OR derivative. Thus d_alpha Psi
vanishes for ANY final OR projection. The final OR calibration uses the true
initial weight, so its normal equations also make d_gamma Psi vanish.

The ordinary initial OR projection and the final calibrated OR projection
need not agree when the outcome model is misspecified. Requiring them to have
the same limit unnecessarily excludes a valid weight-correct branch.

The bounded C2 population audit demonstrates this directly. The target-arm
ordinary OR is projected under the treatment-conditional target distribution;
the weight-correct source calibration projects under the target covariate law.
Their maximum coefficient differences are 0.024343 (control) and 0.025468
(treated). Nevertheless the final score derivatives are at most 5.11e-15,
and finite-difference derivative errors are below 5.56e-12. The identity holds
at both source shifts in the default DGP and agrees under 16/24-point
quadrature. The point bias is zero to the reported numerical precision.

## Correct outcome model

For a transportable source, the initial and final OR limits can both equal
the true conditional outcome model. Then d_gamma Psi vanishes for any final
weight limit. The calibrated weight equations use the true initial OR
derivative, which makes d_alpha Psi vanish as well. Initial and calibrated
weight projections need not be globally identical in this branch.

Both arguments rely on the relevant true model being represented and on the
calibrated losses having the stated population minimizers. Common feature
spaces are needed: OR calibration equations must span the weight-score
derivative directions, and conversely for weight calibration. This is why
the former mismatched-basis simulation could not simply inherit the proof.

## Target anchor must match the implemented method

For the Hou-style target anchor, write the inverse propensity as 1+o_gamma,
where o_gamma is the odds. The OR calibration weights are A*o_init. If the
propensity model is correct, both odds limits are the true odds, and the
calibrated OR equations E[(1-pi) phi (m_true-m_cal)]=0 make the propensity
derivative of the target score vanish. The OR derivative vanishes because
E[A/pi|X]=1. If the OR is correct, the outcome residual makes the propensity
derivative zero, while derivative-weighted PS calibration makes the OR
derivative zero. Again the initial and final misspecified projections need
not both agree.

This argument applies to the calibrated target anchor used in the current
experiments. It does not automatically justify the ordinary AIPW comparison
anchor. The reviewed JASA version still described ordinary target AIPW with
both-nuisance consistency assumptions. Its method description and target
proof should be aligned with the actual Hou-calibrated study before claiming
coverage of the target-misspecified C2/C3 branches.

## What the manuscript proof still needs

Replace blanket agreement of all initial/final pseudo-targets by explicit
branch-specific consistency and calibration identities. Then establish the
appropriate nuisance rates and empirical-process bounds. The remainder in a
misspecified branch can include a squared error of the correctly specified
nuisance as well as a cross-product; a product-only bound is not automatic.
The two-layer training-moment argument likewise needs those actual rates.

This is a population orthogonality check, not a finite-sample inference theorem.
It does not remove the biased-source issue found in the next-paper audit:
with a shifted OR-correct source, one-round target initialization may differ
from that source's outcome truth. Under weight misspecification, regularity
around such a biased candidate still needs separate justification.

Evidence: scratch implementation/next_paper/v9/union_orthogonality_audit/;
script diagnosis/next_paper/probe_union_orthogonality.R. No estimator, DGP or
remote manuscript was changed for this audit.

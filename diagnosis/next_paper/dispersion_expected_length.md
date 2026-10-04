# Coverage-preserving fallback and expected length

The Gaussian adaptation rate is not asserted as new; see
`adaptive_gaussian_literature_update.md` for the relevant 2026 result.

The frozen v8 reference always returned the known-valid reference interval if
its intersection with the dispersion interval was empty. This preserves
coverage, but is unsuitable for claiming an all-valid expected-length rate:
the reference can miss the truth with fixed positive probability, and a
fallback of length O(n^(-1/2)) on that event can dominate a shrinking-in-K
dispersion interval.

Let I_D and I_R have failure budgets alpha_D and alpha_R. Return their
intersection when nonempty, and the SHORTER of them when empty. The event
that both contain the truth implies a nonempty intersection, hence coverage
is still at least 1-alpha_D-alpha_R. For every realized dataset, output length
is at most min(length(I_D),length(I_R)). This simultaneously retains the
reference's global length bound and the dispersion interval's local bound.
The choice is implemented explicitly; the old reference fallback remains an
option. Current-paper RoCE inference and frozen v8 results are unchanged.

## A Gaussian all-valid expected-length upper bound

At the all-valid law, Q has a central chi-square distribution with d residual
degrees of freedom. Write t=log(1/delta). A noncentral chi-square lower-tail
bound follows directly from its negative log-mgf:

    P_lambda[Q <= d+lambda-2 sqrt((d+2lambda)t)] <= exp(-t).

For x=(Q-d)_+, choose L=2x+4 sqrt(dt)+8t. The inequality
2 sqrt(2tL) <= L/2+4t implies
L-2 sqrt((d+2L)t) >= x. Thus the exact CDF-inverted upper bound U_delta(Q)
is no greater than L. Since E(Q-d)_+ <= sqrt(2d),

    E U_delta(Q) <= 2 sqrt(2d)+4 sqrt(dt)+8t,
    E sqrt(U_delta(Q)) <= sqrt(E U_delta(Q)) = O(d^(1/4)).

With fixed validity fractions, bounded homogeneous/private variances of order
1/n, and no shared TATE variance floor, the pooled variance is O(1/(nK)) and
the worst-subset factor C is O(1/(nK)). Since d=O(K),

    E length(I_D) = O(1/sqrt(nK) + 1/(sqrt(n) K^(1/4))).

The shorter-fallback rule inherits this upper bound while its output never
exceeds the separately valid reference length. In the independent Gaussian
submodel this matches the order of the all-valid adaptation lower bound in
majority_oracle_adaptation_bound.md. It is a pointwise adaptation statement
under uniform coverage, not an oracle-efficiency claim over all heterogeneous
configurations. A nonzero shared-target TATE variance adds its own floor.

The result applies to known Gaussian covariance. Estimated covariance,
high-dimensional nuisance fits and an observed-data upper bound still require
their own arguments. A short interval after a wrong validity-count assumption
does not retain this coverage guarantee.

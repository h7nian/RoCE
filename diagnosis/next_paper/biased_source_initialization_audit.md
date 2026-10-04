# Initial OR models and regularity of biased source candidates

This audit checks the population score at the source-specific transported
mean mu_(t,s), which may differ from the target mean mu_t. It does not call
an invalid source valid for the target estimand.

For binary outcomes and a fixed merged weight q, the OR-parameter derivative is

    E_t[h_source(X)*phi(X)] - E_s[I(A=a)*q(X)*h_source(X)*phi(X)],
    h_source(X)=m_source(X)*(1-m_source(X)).

The final weight calibration instead balances the feature vector multiplied
by h_initial. These derivatives agree at the population limit if the initial
OR is correct for that source, or if the final merged weight is truly correct.

In a valid source with the target outcome regression, target initialization
can satisfy the OR condition. At a shifted source in C3, the OR model remains
correct but the merged-weight model is misspecified. The initial target OR
then differs from the source OR. Solving the h_target calibration equation
does not generally solve the h_source equation needed for the biased mean's
orthogonality.

## Population check using the frozen bounded DGP

`probe_biased_source_orthogonality.R` uses the existing four active slopes,
two source positions, radius12, C1/C3 and treated outcome shifts0/.5/1/2.
It solves unpenalized population weight calibration equations and verifies
the resulting OR derivatives by finite differences. The16/24 quadrature
comparison differs by at most5.31e-15.

| C3 source | Fixed shift | Target-initialized score derivative, maximum absolute coordinate |
|---:|---:|---:|
| 1 | .5 | .001999 |
| 1 | 1 | .002864 |
| 1 | 2 | .002273 |
| 2 | .5 | .003007 |
| 2 | 1 | .004315 |
| 2 | 2 | .003434 |

The target-initialized calibration equation itself has residual below1e-14.
Thus the nonzero score derivative is not an observed optimizer failure.
Source-initialized versions have score derivatives below7e-15. Both versions
still reproduce the correct source-specific transported mean at the true
source OR. This is a regularity/orthogonality issue, not nonidentification
of that mean. C1 and zero-shift C3 controls have derivatives at numerical zero.

## Consequences and limits

The current paper does not automatically become invalid: its valid-source
theorem does not assert Gaussian inference for every incompatible candidate.
The proposed dispersion pivot, however, needs appropriate distributional
control for all participating candidates, including biased ones. That stronger
requirement needs source-local initialization, an augmented influence expansion,
or another justified construction in the nonlinear-OR/misspecified-weight branch.

Two-round source-local initialization removes the demonstrated population
derivative mismatch. It does not itself prove all remaining uniform rates,
estimated covariance accuracy or growing-K Gaussian approximation conditions.

The audit uses fixed outcome shifts. Along a genuine local sequence rho_n->0,
the mismatch can vanish. A schematic remainder contribution is of order
rho_n times the OR estimation error; for rho_n=O(n^(-1/2)), this can be negligible
under suitable uniform rates. Therefore this audit does not prove that one-round
inference is impossible for every local-departure model. Fixed rho=.5 must not
be reinterpreted as such a sequence.

For an identity-link linear OR, h is constant, so this particular initialization
mismatch is absent. For a correct merged-weight model, balancing holds for any
integrable h and the mismatch is also absent. Keep these branches distinct.

The default DGP, current-paper estimator and existing campaign configurations
were not changed by this audit. Results are under
`implementation/next_paper/v7/biased_source_population_audit_v1/` on scratch.

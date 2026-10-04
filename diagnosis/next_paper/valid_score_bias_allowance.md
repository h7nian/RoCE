# A Gaussian interval extension allowing valid-score nuisance bias

Sufficient Gaussian construction with an opt-in research implementation.
It addresses the fitted mean drift seen in the90-repeat diagnostic. Covariance
is assumed known, and valid-coordinate mean-error envelopes must be justified
externally. Neither assumption is supplied by this note for the fitted method.
No novelty claim is made.
If the training-data envelopes hold only with high probability, their failure
probability needs its own allocation in the total interval error budget.

Write E[Z]=H theta+r and Var(Z)=Sigma, with two arm means and TATE contrast
ell=(1,-1). A truly valid coordinate set V contains the guaranteed-valid set A
(at least the two target anchors). Allow |r_i|<=b_i for valid coordinates;
invalid-coordinate errors remain unrestricted. Bounds assigned to a candidate
coordinate only apply when that coordinate is valid.

Let c_all be full GLS TATE coefficients and v_all their variance. Let R_A be
the residual precision for the known reference, embedded into the full space:

    R_A = embed[ Sigma_A^-1
           - Sigma_A^-1 H_A (H_A' Sigma_A^-1 H_A)^-1 H_A' Sigma_A^-1 ].

It satisfies R_A H=0 and R_A Sigma R_A=R_A. Define P=I-R_A Sigma. For any
comparison coefficients d supported on V with H'd=ell, the projected
coefficients d_tilde=P d remain supported on V, preserve H'd_tilde=ell and
have variance no larger than d'Sigma d. They are covariance-orthogonal to
the reference residuals.

The difference d_tilde-c_all consequently belongs to the row space of the
existing reduced bias map M. Here M H=0 and M Sigma R_A=0; both properties
were checked for the implemented map. Its uncertain-coordinate block is the
identity, so rank(M) equals the total dimension minus |A|. The constraints
H'x=0 and R_A Sigma x=0 have that same null-space dimension; their ranks add
because H and Sigma times the range of R_A are independent. This identifies
the row space. Also Sigma c_all=H(H'Sigma^-1 H)^-1 ell, so every unbiased
comparison has covariance v_all with full GLS. In the
known-covariance Gaussian model,

    Q=(M Z)' (M Sigma M')^-1 (M Z)

is noncentral chi-square, with Lambda=(M r)'(M Sigma M')^-1(M r), even when
the known reference has nuisance drift. Cauchy–Schwarz then gives

    |c_all' r| <= |d_tilde' r|
                 + sqrt((d_tilde'Sigma d_tilde-v_all) Lambda).

This corrects the missing common-mean allowance; it does not assume that
residual dispersion can detect a common shift.

## Two usable comparison bounds

Reference GLS gives the pair

    B_ref = sum_i |c_ref,i| b_i,
    C_ref = v_ref-v_all.

For shared-prediction covariance, inflate private2x2 blocks as in
shared_covariance_bound.md and take d as the GLS coefficients for an unknown
valid subset under the inflated covariance. Then

    d_tilde'Sigma d_tilde <= d'Sigma d <= v_V(Sigma_bar).

Because P only adds coordinates from A, the transformed envelopes
b_tilde=|P|'b are valid for every coefficient actually used by d. Therefore
|d_tilde'r|<=sum_j |d_j| b_tilde_j, without restricting invalid means.

The shared-prediction GLS coefficients reduce to two anchor coefficients
a=ell-beta and two source-pool coefficients beta. Pool weights within each
arm are positive inverse-private-variance weights. With target covariance
blocks A0 (anchors), B0 (anchor/common), and C0 (common predictions), and
source-pool private variances u1,u0,

    beta = [A0+C0-B0-B0' + diag(u1,u0)]^-1 (A0-B0') ell.

A sufficient envelope is

    B(u1,u0) = sum_a [ |a_a| b_tilde,anchor,a
                   + |beta_a| max_j b_tilde,source_j,a ].

For prescribed positive valid counts, each pool variance lies between the
values obtained from the smallest and largest allowed private variances.
For fixed u0, all coefficients above are affine functions of the same
linear-fractional transform of u1, with positive denominator. A weighted
absolute-value sum is convex along that line segment, so its maximum occurs
at an endpoint. Repeating for u0 bounds B by its four rectangle corners.
This needs positive private variances and a positive-definite displayed
matrix throughout the rectangle. An arm with zero guaranteed sources fixes
its pool coefficient at zero and needs a separate one-dimensional calculation.

Together with the existing worst-subset inflated variance, this supplies

    B_shared = maximum corner envelope,
    C_shared = max_V v_V(Sigma_bar)-v_all.

Keep each B paired with its own C. Choosing the smaller variance bound and
an unrelated smaller bias envelope is not justified. For any upper confidence
bound Lambda_U, the valid deterministic bias allowance is

    min( B_ref + sqrt(C_ref Lambda_U),
         B_shared + sqrt(C_shared Lambda_U) ).

Add that allowance to the full-GLS noise radius. The reference interval used
for intersection must itself include B_ref. All-known-valid and anchor-only
branches likewise require their own GLS bias allowance. At b=0 the construction
should recover the existing bounds and interval behavior exactly.

## Implementation and remaining work

`make_joint_dispersion_calibration()` now accepts `valid_score_bias`: NULL or
zero preserves the previous program exactly; a nonnegative scalar applies
uniformly, and a full vector or (target+sources)-by-two matrix supplies separate
coordinate bounds. Matrix columns, when named, must be mu1/mu0 and rows follow
the covariance's target/source order. The bounds concern conditional score
means, not outcome-model prediction errors or source departure parameters.

All four existing bound modes are supported. Exact subset search stores each
paired bias/variance bound and takes their maximum before comparison with the
reference. Shared-prediction modes use the four-corner envelope and its
uncapped inflated-subset variance. The reference-only mode uses its own pair.
The shared envelope currently requires positive private variances in active
arms; the reference-only mode remains available outside that condition.

835 assertions passed, including the original tests and actual fitted
covariance fixtures. Independent geometry checks enumerated6696 subsets in252
covariance/count cases and exposed531 failures of the unprojected Cauchy step.
API corner bounds matched the independent calculation;16 cases matched the
previous frozen zero-budget calibration and interval outputs exactly. The
standalone tests explicitly compare numeric identities without irrelevant
names; an earlier name-attribute failure is retained in the scratch audit.

Population-quadrature
errors can provide oracle diagnostic budgets, but they are unavailable to a
feasible estimator. Uniform nuisance rates or valid data-dependent envelopes,
estimated-covariance uncertainty and a suitable score-distribution approximation
remain necessary. This extension is a proposed component of that proof route,
not a completed confidence procedure for high-dimensional RoCE. Evidence is
under scratch implementation/next_paper/v10/valid_score_bias_budget_v2/.

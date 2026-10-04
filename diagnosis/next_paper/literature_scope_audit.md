# Scope audit: Guo's papers and the proposed target-anchored extension

Checked against primary publication records on 2026-09-21. This is a research
comparison, not a manuscript draft or novelty claim.

| Reference | Result relevant here | Boundary to retain |
|---|---|---|
| Guo et al., JRSSB 2018, [DOI](https://doi.org/10.1111/rssb.12275) | TSHT with voting and a plurality identification rule | Selection/oracle conclusions require their model and separation conditions |
| Guo, JRSSB 2023, [published article](https://academic.oup.com/jrsssb/article/85/3/959/7174915) | Searching/sampling accounts for validity-selection errors | Condition 3 concerns strongly relevant valid IVs; its finite-sample plurality alternative has further requirements |
| Guo et al., JASA 2025, [published article](https://www.tandfonline.com/doi/abs/10.1080/01621459.2024.2443246) | RIFL provides coverage despite mistakes in selecting the prevailing group | Assumption 1 requires a majority; Remark 6 fixes the site count in the proved asymptotic theory |

The published RIFL Remark 6 distinguishes an algorithm that permits increasing
site count from a proved growing-site asymptotic result. Equations (21)–(23)
give sampling-error exponents 1/[d L(L-1)/2] and 1/[L(L-1)/2]. These expose
potential resampling cost as L grows, but do not establish a growing-L theorem
or an unavoidable complexity lower bound. Theorem 2 gives coverage in its
stated iterated n,M limits. Theorem 3's oracle comparison adds separation and
a condition on the prevailing group's size. These statements replace reliance
on the earlier arXiv version alone.

## Implications for our work (our interpretation)

The designated target already identifies TATE. A guaranteed valid-source count
is an extra condition for borrowing, not the sole identification mechanism.
Source-assisted candidates share target observations; independent-site RIFL
variance formulas cannot be transplanted without accounting for this dependence.
Likewise, independent pairwise perturbations can violate the exact dependence
of differences constructed from common candidate estimates.

Our present reference explicitly handles a shared target component and need
not classify valid sources. This is not enough for a contribution: RIFL already
addresses selection-robust federated inference, and Guo 2023 already uses test
inversion. Necessary further work includes uniform nuisance expansions,
growing-K conditions, estimated covariance, communication requirements and
interval length under stated valid-count assumptions.

The new half-valid Gaussian lower bound clarifies one limit on oracle-relative
length; it is not asserted as a result from these papers or as a full RoCE
lower bound. Broad claims that weak deviations are simply “unidentifiable”
should be replaced by the precise selection, coverage or efficiency statement.

Publisher PDF requests returned HTTP 403 in this environment. The RIFL
published statements above were checked in the publisher-indexed article text,
including numbered assumptions, theorems and remarks; the older author PDF
remains separately labeled as a preprint in scratch. No inaccessible appendix
has been represented as newly verified.

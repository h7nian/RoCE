# Relevant 2026 result and the remaining research scope

Primary source checked on 2026-09-21:
[Wang, Chai and Gao, Adaptive Confidence Intervals in Efron's Gaussian Two-Groups Model](https://arxiv.org/abs/2604.26992).

They characterize adaptive interval length for independent Gaussian mean-shift
contamination with unknown contamination proportion. With known common noise
variance, the all-valid adaptation term is sigma K^(-1/4); their full rate also
depends on the actual contamination proportion. They give a polynomial-time
Fourier certification construction. Unknown variance creates further adaptation
costs. This is a preprint, and its model and guarantees should be cited rather
than conflated with our setting.

## Consequences for our work (our assessment)

The Gaussian K^(-1/4) obstruction is not a new standalone contribution. Our
earlier calculations remain useful checks and extensions with a trusted anchor,
but must not be marketed as establishing a previously unknown Gaussian rate.
We also have not proved globally optimal length over all biased configurations;
our target-capped dispersion intervals can remain target-scale away from the
all-valid region.

The relevant directions to establish are:

- within-site replication and observed-data causal likelihoods;
- unknown, heterogeneous site noise, with local variance information;
- shared target contributions and their variance floor;
- different valid source sets in the two treatment arms;
- high-dimensional calibrated nuisance fits and communication constraints;
- minority-valid borrowing made possible by a trusted target, with its limits.

These are candidate distinctions requiring proofs and comparisons, not novelty
claims. The bounded replicated-score construction and binary likelihood lower
bound address simpler observed-data submodels. Their extension to fitted
high-dimensional source-assisted scores remains open.

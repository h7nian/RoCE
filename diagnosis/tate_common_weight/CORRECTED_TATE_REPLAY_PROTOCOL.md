# Full original-estimator replay after the CV normalization correction

Use the verified corrected installed package, saved seed10013 C1/rho0 data
and the exact original five outer-fold partitions. Call `run_tate_crossfit`
directly, with one-round communication, both arms, both sources, 100 nuisance
lambdas, min selection, aggregation cutoff1, and M_fit=M_inference=5.
No experimental projection or shared-final-basis modification is included.

Arms and sources run sequentially so R warning conditions are retained at the
parent boundary; native nuisance CV uses five allocated threads. The direct
entry point has no bootstrap argument or comparison-method bootstrap path.
Preserve complete fitted TATE state, warnings and input provenance. Verify
point reconstruction from TATE pseudovalues, site-centered variance
reconstruction, and unchanged target-only estimate against the old fit.

This is a fixed-data replay, not a new independent replication or a coverage
gate. It isolates the implemented CV normalization change within the full
original estimator. Full warning provenance is still explicitly incomplete;
the replay's SE must not be presented as a validated final inference method.

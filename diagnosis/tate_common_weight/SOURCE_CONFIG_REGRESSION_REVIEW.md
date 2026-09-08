# Configuration-aware source selector: regression result

The first C1 regression job, 18444629, failed at the end-of-run code-hash
guard because its included protocol document was edited while it ran. No
result bundle was committed. Preserve its Slurm logs; this is an input
stability failure, not evidence of a numerical solver failure.

Replacement job 18444788 completed with exit code 0 in 50 seconds after
the referenced files were left unchanged. In `source_stagewise_selection_v2/`,
the candidate loss CSV, complete selection-state RDS and outer-fold report
CSV are each byte-identical to v1. All five v2 payload checksums passed.
Only configuration/provenance metadata differ as intended.

C2/C3 selection array 18444632 additionally depends on this successful
regression, remaining validation audits 18443708, and full outer-training
audits 18444628. These successful C1 checks do not certify uncompleted C2/C3
jobs or the full TATE method. Leave files listed in the in-flight jobs'
code-hash manifests unchanged; record later status in separate review notes.

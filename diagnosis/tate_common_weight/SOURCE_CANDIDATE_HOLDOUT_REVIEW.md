# First-split source candidate evaluation

The training-noise coupled candidates were evaluated on global validation
fold 2, using the audited outer-1-excluded nuisance fits. Both arms and all
three scales 0.5,1,2 completed without candidate failure. Corrected means
reconstructed as score minus the full moment/coefficient product to 1e-10.

| Arm | Original mean | Corrected, scale 0.5 | Corrected, scale 1 | Corrected, scale 2 |
| --- | ---: | ---: | ---: | ---: |
| Control | 0.5292658 | 0.5166152 | 0.5263259 | 0.5292658 |
| Treated | 0.7309833 | 0.6853606 | 0.7239084 | 0.7309833 |

The bundle `source_validation_candidate_1_2_v1/` under the independent-pilot
root retains full paired site observation vectors, all six candidate fits,
48 conditional block risks, and original IDs. All four payload checksums
passed. No nuisance refit or scale selection occurred.

Initial-block risks use different downstream coefficient values across these
scale candidates. They are explicitly labeled as such and must not be
summed into a purported joint convex risk or compared as though downstream
parameters were fixed. Final-block selection and conditional initial-block
selection need a consistently enumerated candidate policy and all validation
folds. These one-fold means cannot justify selecting a scale or reporting a
final TATE interval.

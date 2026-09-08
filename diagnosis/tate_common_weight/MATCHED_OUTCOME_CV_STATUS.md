# Matched high-dimensional CV task status

Array 18453175 failed both tasks before fitting because the new driver had
an extra closing parenthesis in its warning-capture wrapper. No CV/input
bundle was committed. Both original Slurm error logs remain preserved.

The driver syntax was corrected and `parse()` now passes. Both intended
output directories were confirmed absent before resubmission. Data, explicit
fold seed, lambda grid, nuisance weights, package versions and CPU allocation
are unchanged. A resubmitted task is not evidence of numerical success;
check its result and input identity before comparing CV paths.

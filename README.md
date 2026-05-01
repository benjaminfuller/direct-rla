# direct-rla
Direct Ballot Selection Risk Limiting Audits


This repository generates samples size calculations for ballot comparison audits, polling using Minerva, and direct ballot selection.

For a mode in direct, comparison, polling you want to run the following ./dup_for_diff_sizes.sh mode.  The options for mode are polling, comparison, and direct. Polling and
comparison take a couple of minutes to run. Direct takes a couple of days to run. As desired one can edit the script to reduce the number of sizes being run or the sizes. For population N, time scales
with Nlog N.

Under the hood this runs three python scripts for the six sizes present in the paper.

1. python3 duplicate_audit_math.py number mode
2. python3 combine_results.py mode
3. python3 process_duplicate_results.py mode

This results in three main files:

1. results_combined_mode.csv the believed best solutions across a variety of parameters
2. time_results_by_sizemode.csv the time savings across sizes for a statistically accurate manifest in comparison to the same method with a full manifest.
3. time_results_by_sizemode.png a graph of the time savings across sizes for a statistically accurate manifest in comparison to the same method with a full manifest.


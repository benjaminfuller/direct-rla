# direct-rla
Direct Ballot Selection Risk Limiting Audits


This repository generates samples size calculations for ballot comparison audits, polling using Minerva, and direct ballot selection.

For a mode in direct, comparison, polling you want to run the following three commands


1. ./dup_for_diff_sizes.sh <mode>
2. python3 combine_results.py <mode>
3. python3 process_duplicate_results.py <mode>

This results in three main files:

1. results_combined_<mode>.csv the believed best solutions across a variety of parameters
2. time_results_by_size<mode>.csv the time savings across sizes for a statistically accurate manifest in comparison to the same method with a full manifest.
3. time_results_by_size<mode>.png a graph of the time savings across sizes for a statistically accurate manifest in comparison to the same method with a full manifest.


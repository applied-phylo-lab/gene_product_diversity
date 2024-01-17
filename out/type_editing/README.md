This directory is for RNA editing type gene product diversity.

out_basic.txt contains steady-state predictions for different parameter combinations.

n_all.txt contains combinations of 3 types of modification events with which genome-wide distributions are simulated. Columns correspond to deleterious, beneficial, and neutral modifications respectively.

Files with "distr_mix_" contain steady-state predictions for modification events with parameter values drawn from pre-specified distributions. Number at the end of each file name corresponds to the number of different types of modifications (as ordered in n_all.txt).

Files with "moment_out_" contain moments of genome-wide distributions. Number at the end of each file name denote product abundance cutoffs (if the modified isoform's abundance is below the cutoff, the modification event is not counted).

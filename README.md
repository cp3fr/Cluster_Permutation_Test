# Cluster_Permutation_Test

Custom Matlab Scripts for Running the Cluster Permutation Test (FieldTrip) on EEG

How to use the script:
- download fieldtrip from http://www.fieldtriptoolbox.org
- download eeglab from https://sccn.ucsd.edu/eeglab/download.php
- place some sample EEG data into the same folder where the scripts are
- in batch.m edit the variables: 'path_to_eeglab', 'path_to_fieldtrip', and 'pathname_inputfile'
- run the script
- it should give as ouput figures showing the clusters and their statistics

Note: Currently the script runs an analysis comparing single-trial EEG for from two conditions of one subject.
In order to run across-subjects analysis the scripts need modification. That is, compute conditon-average ERPs 
for each subject, and add those data to 'data' variable. 
In 'run_cpt.m' line 128, cfg.satistics must be set to 'ft_statfun_depsamplesT', because the samples are dependent.

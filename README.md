# Cluster_Permutation_Test

Custom Matlab Scripts for Running the Cluster Permutation Test (FieldTrip) on EEG

Note: Currently the script runs an analysis comparing single-trial EEG for from two conditions of one subject.
However, input data should be conditon-average ERPs from two conditions from multiple subjects.
(needs adaptation)

How to use the script:
- download fieldtrip from http://www.fieldtriptoolbox.org
- download eeglab
- place some sample EEG data into the same folder where the scripts are
- in batch.m edit the variables: 'path_to_eeglab', 'path_to_fieldtrip', and 'pathname_inputfile'
- run the script
- it should give as ouput figures showing the clusters and their statistics

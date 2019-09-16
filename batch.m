%%batch.m
clear all;close all;clc;

%path settings
path_to_eeglab = '/Volumes/methlab/NLP/Ce_ETH/Christian/UnfoldNeuroimage/toolboxes/eeglab/';
path_to_fieldtrip = '/Volumes/methlab/NLP/Ce_ETH/Christian/UnfoldNeuroimage/toolboxes/fieldtrip-20190218/';
pathname_inputfile = '/Users/chpfei/Repos/Cluster_Permutation_Test/sample_data/Struct.mat';

%add toolboxes
addpath(path_to_eeglab); eeglab; close;
addpath(path_to_fieldtrip); ft_defaults;

%load alldata
load(pathname_inputfile,'EEG');
alldata = EEG;
clear EEG;

%which conditions to compare
conditions_to_compare = [2 1]; %'Sync', 'Async'

%sampling rate
sampling_rate = alldata(1).srate;

%eeg channels only
chans = 1:19;

%channel location file
chanlocs = alldata(1).chanlocs;

%times
times = alldata(1).times;

%number of permuations
nperm = 100; %1000 is better

%% Note: Here the analysis is run on single trials for one subject
%        Typically you would like to run it conditon-average ERPs from multiple subjects

%prepare the data
data = cell(1,length(conditions_to_compare));

%loop over conditions
for i = 1:length(conditions_to_compare)

  %extract the EEG data only (drop heart channel)
  eeg_data = alldata(conditions_to_compare(i)).data(chans,:,:);

  %cell array of conditions, containing double array of channel x time x trials
  data(1,i) = {eeg_data};

end



%run the cluster permuation test
clust = run_cpt(data,... %chan x sp x cond x subj
  struct(...
    'srate',sampling_rate,...
    'chanlocs',chanlocs,...
    'times',times,...%by default the test is ran across all sampling points, but could be reduced to samples of interest or averages across windows
    'nperm',nperm,...
    'alpha',0.05)); %chan x sp x cond x subj  st.(mn).epoDiff.ci =computeCiWithin(permute(shiftdim(st.(mn).epoDiff.data,-1),[2,3,1,4]),struct('alpha',st.config.alpha,'dimLoop',[1 2])); %chan x sp x cond x ci



%plot the results
for vn = {'pos','neg'}

  if ~isempty(clust.(vn{1}))

    figure
    imagesc(clust.(vn{1}).chanTime)
    title(sprintf('%s cluster stats: t=%f, p=%f',vn{1},clust.pos.t,clust.pos.p))
    xlabel('sampling points')
    ylabel('channels')

  end

end
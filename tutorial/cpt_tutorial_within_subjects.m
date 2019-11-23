%cpt_tutorial_within_subjects.m
%
%based on:
%http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/

clear all;close all;clc;
restoredefaultpath;

%add toolbox
addpath('/Users/chpfei/Repos/fieldtrip-20190218/'); 
ft_defaults;

%load data
load ERF_orig;


cfg = [];
cfg.channel = {'MEG'};
cfg.latency = [0 1];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;

%neighbors structure
cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, allsubjFC{1,1});
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;

%use the latest spm version, because the default 'spm8' causes a mex file error                  
%see here for more details about the workaround:
%http://www.fieldtriptoolbox.org/faq/matlab_complains_about_a_missing_or_invalid_mex_file_what_should_i_do/
cfg.spmversion     ='spm12';

%design structure
subj = 10;
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
cfg.design = design;

cfg.uvar  = 1;
cfg.ivar  = 2;

%run cluster permutation test
[stat] = ft_timelockstatistics(cfg, allsubjFIC{:}, allsubjFC{:})

%save results
save stat_ERF_planar_FICvsFC_GA stat

%load results
load stat_ERF_planar_FICvsFC_GA

% load individual subject data
load('ERF_orig');
% calculate the grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_FC         = ft_timelockgrandaverage(cfg,allsubjFC{:});
GA_FIC        = ft_timelockgrandaverage(cfg,allsubjFIC{:});
% "{:}" means to use data from all elements of the variable

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_FICvsFC = ft_math(cfg,GA_FIC,GA_FC);

figure;
% define parameters for plotting
timestep = 0.05;      %(in seconds)
sampling_rate = allsubjFIC{1,1}.fsample;
sample_count = length(stat.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
% get relevant (significant) values
pos_cluster_pvals = [stat.posclusters(:).prob];

% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(GA_FICvsFC.label, stat.label);

% plot
for k = 1:20;
   subplot(4,5,k);
   cfg = [];
   cfg.xlim=[j(k) j(k+1)];
   cfg.zlim = [-5e-14 5e-14];
   pos_int = zeros(numel(GA_FICvsFC.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout = 'CTF151_helmet.mat';
   ft_topoplotER(cfg, GA_FICvsFC);
end
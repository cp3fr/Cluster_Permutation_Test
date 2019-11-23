clear all;close all;clc;
restoredefaultpath;

%add toolboxes to path
addpath('/Users/chpfei/Repos/fieldtrip-20190218/'); 
ft_defaults;
addpath('/Users/chpfei/Repos/Cluster_Permutation_Test/')

%load input data
load chanlocs;

%load and prepare input data
%[subj x cond] cell array with [chan x sp] double array, condition-average EEG
load data3; 
data={};
for icond = 1:length(data3)
  for isubj = 1:size(data3{icond},3)
    data(icond,isubj)={squeeze(data3{icond}(:,:,isubj))};
  end
end;
clear data3;

%cluster permutation test within-subjects
cfg = [];
cfg.alpha = 0.05;
cfg.nperm = 500;
cfg.srate =  500;
cfg.minchan = 2;
cfg.method = 'distance';
[s] = clusterperm_within(data, chanlocs, cfg)

%inspect results

n = [length(s.pos.p), length(s.neg.p)];
panels = [2,max(n)];
clust = {'pos','neg'};
figure;

for iclust = 1:length(clust)
  count = max(n) * (iclust-1);
  for i = 1:n(iclust)

    count = count+1;

    subplot(panels(1),panels(2),count);

    imagesc(s.(clust{iclust}).chanTime==i)

    title(sprintf('%s cluster\nt = %.4f\np = %.4f\nh = %d', clust{icond}, s.(clust{iclust}).t(i),s.(clust{iclust}).p(i),s.(clust{iclust}).h(i)    ));
  end
end

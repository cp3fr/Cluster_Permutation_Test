function [s] = run_cpt(data,varargin)
%function [s] = run_cpt(data,cfg)
%
% Cluster permutation test for BETWEEN TRIALS DESIGN as implemented in fieldtrip 
%
% data:        1xcond cells of chan x sp x trial EEG data
% cfg:         configuration file (optional) with fields
%    .alpha:   significance threshold (default: 0.05)
% s:           statistical results
%
% based on http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
% Currently only working between two conditions.. 
% Implement analysis across more conditions too


%convert data to fieldtrip structure
%perform cluster permutation test
%find significant clusters
%make output structure with significant electrode count per cluster and
%pvalues


%% Configuration file

if length(varargin)==1
  
  CFG = varargin{1};
  
else
  
  CFG = [];
  
end

if ~isfield(CFG,'alpha')

  CFG.alpha = 0.05;

end

if ~isfield(CFG,'srate')
  
  CFG.srate = 125;
  
end

if ~isfield(CFG,'chanlocs')

  load('eegchanlocs.mat','chanlocs');

  CFG.chanlocs = chanlocs;

  clear chanlocs;

end

CFG.nCond = length(data);

[CFG.nChan,CFG.nSp,CFG.nTrials] = size(data{1});

if ~isfield(CFG,'times')

  CFG.times = 1:CFG.nSp;
  
end

if ~isfield(CFG,'nperm')
  
  CFG.nperm = 1000;
  
end


%% Prepare the data

d=cell(CFG.nCond,1);

for iCond = 1:CFG.nCond
  
  EEG = eeg_emptyset; %make a EEG structure
  
  EEG.data = squeeze(data{iCond}); %chan x sp x subjs
  
  EEG.nbchan = size(EEG.data,1);
  
  EEG.trials = size(EEG.data,3);
  
  EEG.srate = CFG.srate;
  
  EEG.pnts = size(EEG.data,2);
  
  EEG.xmin = CFG.times(1);
  
  EEG.xmax = CFG.times(end);
  
  EEG.times = CFG.times;
  
  EEG.chanlocs = CFG.chanlocs;
  
  ftPreprocess = eeglab2fieldtrip( EEG, 'preprocessing', 'none' ); %convert eeglab EEG structure of the current condition to fieldtrip data

  cfg = [];
  
  cfg.keeptrials = 'yes';
  
  ftTimelock  = ft_timelockanalysis(cfg, ftPreprocess); %convert to FT timelock structure

  d(iCond) = {ftTimelock};
  
  clear EEG ftPreprocess ftTimelock cfg;
  
end


%% Clusterpermutation test

%configuration structure
cfg = [];

% use the Monte Carlo Method to calculate the significance probability
cfg.method = 'montecarlo';       

% use the independent samples T-statistic as a measure to evaluate the 
%effect at the sample level
cfg.statistic = 'ft_statfun_indepsamplesT'; 

% how to correct..
cfg.correctm = 'cluster';%'no','fdr','bonferoni','max','holms'

% alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clusteralpha = CFG.alpha;    

% test statistic that will be evaluated under the permutation distribution.
cfg.clusterstatistic = 'maxsum'; 

% minimum number of neighborhood channels that is required for a selected 
% sample to be included in the clustering algorithm (default=0).
cfg.minnbchan = 2;           

% see below
% cfg.neighbours = neighbours;   % see below

% -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.tail = 0;

cfg.clustertail = 0;

% alpha level of the permutation test
cfg.alpha = CFG.alpha / 2;         

% number of draws from the permutation distribution
cfg.numrandomization = CFG.nperm;      

%design vector
design = [];

%loop over conditions
for iCond = 1:CFG.nCond

  %get number of events for current condition
  nEvents = size(d{iCond}.trial,1);

  %add condition labels to the design vector
  design(1,end+1:end+nEvents) = ones(1,nEvents).*iCond;

end

% design matrix
cfg.design = design;        

% number or list with indices indicating the independent variable(s)
cfg.ivar  = 1;                   

%make the channel neighbor structure
cfg_neighb        = [];
cfg_neighb.method = 'triangulation';   
neighbours        = ft_prepare_neighbours(cfg_neighb, d{1});

% the neighbours specify for each sensor with which other sensors it can form clusters
cfg.neighbours    = neighbours;  

% cell-array with selected channel labels
cfg.channel       = {'EEG'};    

% time interval over which the experimental conditions must be compared (in seconds)
cfg.latency       = [d{1}.time(1), d{1}.time(end)];       

%use the latest spm version, because the default 'spm8' causes a mex file error                  
%see here for more details about the workaround:
%http://www.fieldtriptoolbox.org/faq/matlab_complains_about_a_missing_or_invalid_mex_file_what_should_i_do/
cfg.spmversion     ='spm12';

%run cluster permutation test
%..for two conditions
if CFG.nCond == 2
  
  [stat] = ft_timelockstatistics(cfg, d{1}, d{2});
  
%..for three conditions
elseif CFG.nCond == 3
  
  [stat] = ft_timelockstatistics(cfg, d{1}, d{2}, d{3});
  
%..or raise an error
else
  
  error('Error: less than 2 or more than 3 conditions.')
  
end

%% output structure 

s=[];

s.pos = []; %positive cluster statistics

if ~isempty(stat.posclusters)

  s.pos.t_dimord = 'clust';
  
  s.pos.t = [stat.posclusters(:).clusterstat]';
  
  s.pos.p_dimord = 'clust';
  
  s.pos.p  = [stat.posclusters(:).prob]';
  
  s.pos.h_dimord = 'clust';
  
  s.pos.h = s.pos.p < CFG.alpha;
  
  s.pos.chanTime_dimord = 'chan_sp';
  
  s.pos.chanTime = stat.posclusterslabelmat;
  
  s.pos.alpha = CFG.alpha;
  
  s.pos.nperm = CFG.nperm;

end

s.neg = []; %negative cluster statistics

if ~isempty(stat.negclusters)
  
  s.neg.t_dimord = 'clust';
  
  s.neg.t = [stat.negclusters(:).clusterstat]';
  
  s.neg.p_dimord = 'clust';
  
  s.neg.p  = [stat.negclusters(:).prob]';
  
  s.neg.h_dimord = 'clust';
  
  s.neg.h = s.neg.p < CFG.alpha;
  
  s.neg.chanTime_dimord = 'chan_sp';
  
  s.neg.chanTime = stat.negclusterslabelmat;
  
  s.neg.alpha = CFG.alpha;
  
  s.neg.nperm = CFG.nperm;

end

s.stat      = stat;


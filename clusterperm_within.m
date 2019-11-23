%%function [s] = clusterperm_within(data, chanlocs, cfg)
%
% Cluster Permutation Test WITHIN SUBJECTS
%
% data:        ncond x nsubj cell of chan x sp double array 
%              of condition-average EEG
% chanlocs:    EEGlab chanloc structure
% cfg:         configuration file (optional) with fields
%    .alpha:   significance threshold (default: 0.05)
%    .nperm:   number of permutations (default: 500)
%    .srate:   sampling rate of the data (default: 500)
%    .minchan: neighbor minimum number of channels 
%              (default: 2)
%    .method:  neighbor channel method (default: 'distance')
% s:           statistical results
%    .stat:    fieldtrip 'stat' structure
%    .pos:     results for positive clusters (taken from stat)
%        .t:   t-values for all positive clusters
%        .p:   probability values for all positive clusters
%        .h:   whether the cluster shows significant differences
%              between the two conditions
%        .chanTime: channel-time matrix showing cluster 
%                   channels and timing (non-zero values 
%                   indicate the cluster number)
%    .neg:     results for negative clusters (taken from stat)
%              substructures are the same as for s.pos
%
% based on:
% http://www.fieldtriptoolbox.org/faq/how_are_electrodes_magnetometers_or_gradiometers_described/
%
% requires fieldtrip toolbox to be added to the matlab path
%
% Only works on two conditions
%
% christian.pfeiffer@uzh.ch
% 23.11.2019
%
function [s] = clusterperm_within(data, chanlocs, varargin)

  %% GENERAL PROCESSING SETTINGS

  if length(varargin)==1
    CFG = varargin{1};
  else
    CFG = [];
  end

  if ~isfield(CFG,'alpha')
    CFG.alpha = 0.05;
  end

  if ~isfield(CFG,'srate')
    CFG.srate = 500;
  end

  if ~isfield(CFG,'nperm')
    CFG.nperm = 500;
  end

  if ~isfield(CFG,'minchan')
    CFG.minchan = 2;
  end

  if ~isfield(CFG,'method')
    CFG.method = 'distance';
  end

  CFG.times = [1 : size(data{1,1},2)] .* (1 / CFG.srate);


  %% CONVERT DOUBLE ARRAY INPUT DATA TO FIELDTRIP TIMELOCK STRUCTURE

  alldata=cell(size(data));

  for icond = 1:size(data,1)
    for isubj = 1:size(data,2)
      
      timelock = struct();
      timelock.label = {chanlocs.labels}';
      timelock.fsample = CFG.srate;
      timelock.avg = data{icond,isubj};
      timelock.time = CFG.times;
      timelock.dimord = 'chan_time';
      timelock.cfg = struct();

      elec = struct();
      elec.elecpos = zeros(length(chanlocs), 3);
      for ind = 1:length( chanlocs )
          elec.label{ind} = chanlocs(ind).labels;
          if ~isempty(chanlocs(ind).X)
              elec.elecpos(ind,1) = chanlocs(ind).X;
              elec.elecpos(ind,2) = chanlocs(ind).Y;
              elec.elecpos(ind,3) = chanlocs(ind).Z;
          else
              elec.elecpos(ind,:) = [0 0 0];
          end
      end
      elec.pnt = elec.elecpos;
      timelock.elec = elec;

      alldata(icond,isubj) = {timelock};
      
      clear timelock elec ind;
      
    end
  end


  %% CONFIGURATION FILE FOR CLUSTER PERMUATION TEST 

  cfg = [];
  cfg.channel = {'EEG'};
  cfg.latency = [];
  cfg.method = 'montecarlo';
  cfg.statistic = 'depsamplesT';
  cfg.correctm = 'cluster';
  cfg.clusteralpha = CFG.alpha;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = CFG.minchan;

  %neighbors structure
  cfg_neighb        = [];
  cfg_neighb.method = CFG.method;
  neighbours        = ft_prepare_neighbours(cfg_neighb, alldata{1,1});
  cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

  cfg.tail = 0;
  cfg.clustertail = 0;
  cfg.alpha = CFG.alpha / 2;
  cfg.numrandomization = CFG.nperm;

  %use the latest spm version, because the default 'spm8' causes a mex file error                  
  %see here for more details about the workaround:
  %http://www.fieldtriptoolbox.org/faq/matlab_complains_about_a_missing_or_invalid_mex_file_what_should_i_do/
  cfg.spmversion     ='spm12';

  %design structure
  subj = size(alldata,2);
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
  [stat] = ft_timelockstatistics(cfg, alldata{1,:}, alldata{2,:});



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


end









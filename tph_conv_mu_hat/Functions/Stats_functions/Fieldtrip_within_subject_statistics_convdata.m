%% Within group tests FIELDTRIP on convolution data

% Maria:

% (b) See whether we can do within-subject statistics comparing the target
% interval to a baseline period (e.g. [-200, 0] ms) - potentially t-test for
% dependent samples

% If we do (b), then we will average the baseline over 200 values and
% construct a new matrix of identical values (repmat of the avg) and compare
% it to the target interval for within-group stats

% TpH:

% Following the 'within-trials' section of: http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#within-trial-experiments

% BUT NEW PERHAPS COMBINE THIS TUTORIAL WITH THIS ONE: http://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_stats/#1-compute-within-participant-contrasts

% As this one shows how you might be able to grandavaearage first, then select data,
% and then change to run the 'ft_statfun_actvsblT' test?


%% MARIA EMAILS

% SWITCH:
% I only changed case epsi2 in the code, change epsi3 code accordingly
% 
% 
% mycase = case1;
% 
% switch mycase
%     case 'epsilon2 early alpha'
%         frequency = [8 12];
%         sampling_freq = 256; % sampling points
%         baseline_period = 200; % in ms
%         baseline_nsampl = ceil((sampling_freq*baseline_period)/1000);
% 
%         % Baseline
%         for i = 1:21
%             % State anxiety
%             StA_baseline{i,1} = ft_conv_sta_epsi2{i};
%             StA_baseline{i,1}.time  = StA_baseline{i,1}.time([1:baseline_nsampl]);
%             StA_baseline{i,1}.powspctrm = StA_baseline{i,1}.powspctrm(:,:,1:baseline_nsampl);
%             % Controls
%             Cont_baseline{i,1} = ft_conv_cont_epsi2{i};
%             Cont_baseline{i,1}.time = Cont_baseline{i,1}.time([1:baseline_nsampl]);
%             Cont_baseline{i,1}.powspctrm = Cont_baseline{i,1}.powspctrm(:,:,1:baseline_nsampl);
%         end
%         
%         mytime = ft_conv_sta_epsi2{1}.time; %time vector in ms  -200:2000 ms    in one participant
% 
%        
%         TOI = [0.1 0.5]; %in seconds, careful, fieldtrip time field is in seconds
%         activation_period_start = find(min(abs(TOI(1) - mytime)) == abs(TOI(1)  - mytime));
%         activation_period_end = find(min(abs(TOI(2) - mytime)) == abs(TOI(2)  - mytime));
        
%% RUN

rmpath(genpath('E:\tpDATA\MATLAB_Scripts\tanx_osc\SPM\tph_convolution_modelling'));
clear all;
%% Directory
rootFolder = 'E:\tpDATA\MATLAB_Scripts\tanx_osc\SPM\tph_convolution_modelling\tph_convolution_modelling_FOI_8_30b';
cd (rootFolder);
% Add Functions
addpath(genpath([rootFolder '\Functions\']));
addpath(genpath([rootFolder '\Data\']));
% Add SPM
spmpath = 'E:\spm12';
addpath(spmpath);
spm('defaults', 'eeg');
% FIELD TRIP
ftpath = 'E:\fieldtrip-20181024';
% Add FieldTrip
addpath(genpath(ftpath));
ft_defaults;
%% Make a Fieldtrip template for SPM data
load('tanx_wavelet_Cont_win_IQR_AB_baselined.mat')
ft_template = Cont_wavelet_win_AB{1, 1};
rmpath(genpath(ftpath));
%% Load each SPM file and save fttimelock to a FieldTrip template

expgroup = [20 21 23:43];
expgroup = setdiff(expgroup, [36,37]);
contgroup = [1:19 22 37];

condition_win = 1;
condition_lose = 2;
condition_noresp = 3;
condition_epsi2 = 4;
condition_epsi3 = 5;

% expgroup
for i = expgroup
    % Load CONV MODEL continous TF data
    filetag = convertStringsToChars(num2str(i)); % nametag of SPM datafile
    subjfolder = convertStringsToChars(sprintf('participant_%d', i));
    D = spm_eeg_load([rootFolder '\Data\SPM_tf\' subjfolder '\conv_model_atf_d_' filetag '_tosc_spm_cont.mat']);
    % Change Time
    ft_template.time = D.fttimelock.time;
    % Change Freq
    ft_template.freq = D.fttimelock.freq;
    % The DATA
    ft_template.powspctrm = [];
    ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi3,:,:,:));
    ft_conv_sta_epsi3{i,1} = ft_template;
    clear S D
end

% For any missing participants
ft_conv_sta_epsi3 = ft_conv_sta_epsi3(~cellfun('isempty',ft_conv_sta_epsi3));

% contgroup
for i = contgroup
    % Load CONV MODEL continous TF data
    filetag = convertStringsToChars(num2str(i)); % nametag of SPM datafile
    subjfolder = convertStringsToChars(sprintf('participant_%d', i));
    D = spm_eeg_load([rootFolder '\Data\SPM_tf\' subjfolder '\conv_model_atf_d_' filetag '_tosc_spm_cont.mat']);
    % Time
    ft_template.time = D.fttimelock.time;
    % Freq
    ft_template.freq = D.fttimelock.freq;
    % The DATA
    ft_template.powspctrm = [];
    ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi3,:,:,:));
    ft_conv_cont_epsi3{i,1} = ft_template;
    clear S D
end

% For any missing participants
ft_conv_cont_epsi3 = ft_conv_cont_epsi3(~cellfun('isempty',ft_conv_cont_epsi3));

%% Update Paths
rmpath(genpath(spmpath));
addpath(genpath(ftpath));
ft_defaults;

%% BASELINE = New matrix for powspctrm + time

case1 = 'epsilon2 early alpha';
case2 = 'epsilon3 late beta';

clear StA_baseline Cont_baseline StA_activation Cont_activation ...
    grandavg_sta_baseline grandavg_cont_baseline grandavg_cont_activation ...
    grandavg_sta_activation

TOI = case1;

switch TOI
    case 'epsilon2 early alpha'
        frequency = [8 12];
        sampling_freq = 256; % sampling points
        baseline_time = 200; % in ms
        baseline_nsamples = ceil(sampling_freq/1000*baseline_time);
        
        % MARIA
%         baseline_nsampl = ceil((sampling_freq*baseline_period)/1000)

        % Baseline
        for i = 1:21
            % State anxiety
            StA_baseline{i,1} = ft_conv_sta_epsi3{i};
            % Maria - this is sampling points not TIME
%             StA_baseline{i,1}.time  = StA_baseline{i,1}.time([1:baseline_nsampl]);
            StA_baseline{i,1}.time = [1:baseline_nsamples];
            % And here too == StA_baseline{i,1}.powspctrm = StA_baseline{i,1}.powspctrm(:,:,1:baseline_nsampl);
            StA_baseline{i,1}.powspctrm = StA_baseline{i,1}.powspctrm(:,:,StA_baseline{i,1}.time);
            % Controls
            Cont_baseline{i,1} = ft_conv_cont_epsi3{i};
            Cont_baseline{i,1}.time = [1:baseline_nsamples];
            Cont_baseline{i,1}.powspctrm = Cont_baseline{i,1}.powspctrm(:,:,Cont_baseline{i,1}.time);
        end
        
       % TOI == 100-500ms, so that's a baseline = -200-0ms * 2 = 400ms
        toi_start = 100;
        baseline_multiplier = 2;
        activation_period = baseline_nsamples*baseline_multiplier;
        % Activation start = +100 ms
        period_from_start =  200 + toi_start;
        activation_period_start = ceil(sampling_freq/1000*period_from_start);
        activation_period_end = (activation_period_start + activation_period)-1;
        
    case 'epsilon3 late beta'
        frequency = [13 25];
        sampling_freq = 256; % sampling points
        baseline_time = 100; % in ms
        baseline_nsamples = ceil(sampling_freq/1000*baseline_time);
        % Baseline
        for i = 1:21
            % State anxiety
            StA_baseline{i,1} = ft_conv_sta_epsi3{i};
            StA_baseline{i,1}.time = [baseline_nsamples+1:baseline_nsamples*2];
            StA_baseline{i,1}.powspctrm = StA_baseline{i,1}.powspctrm(:,:,StA_baseline{i,1}.time);
            % Controls
            Cont_baseline{i,1} = ft_conv_cont_epsi3{i};
            Cont_baseline{i,1}.time = [baseline_nsamples+1:baseline_nsamples*2];
            Cont_baseline{i,1}.powspctrm = Cont_baseline{i,1}.powspctrm(:,:,Cont_baseline{i,1}.time);
        end
        % TOI == 1400-1700ms, so that's a baseline = -100-0ms * 3 = 300ms
        toi_start = 1400;
        baseline_multiplier = 3;
        activation_period = baseline_nsamples*baseline_multiplier;
        % Activation start = +1400 ms
        activation_period_start = 410; % I need to automate this
        activation_period_end = (activation_period_start + activation_period)-1;
end

for i = 1:21
    % State anxiety
    StA_activation{i,1} = ft_conv_sta_epsi3{i};
    StA_activation{i,1}.time = [activation_period_start:activation_period_end];
    StA_activation{i,1}.powspctrm = StA_activation{i,1}.powspctrm(:,:,StA_activation{i,1}.time);
    % Controls
    Cont_activation{i,1} = ft_conv_cont_epsi3{i};
    Cont_activation{i,1}.time = [activation_period_start:activation_period_end];
    Cont_activation{i,1}.powspctrm = Cont_activation{i,1}.powspctrm(:,:,Cont_activation{i,1}.time);
end


%% STATS using grand avaerage approach

% Grand average = keeping individuals
cfg = [];
cfg.channel = [1:64];
cfg.latency        = ('all');
% For plotting
%    cfg.keepindividual = ('no');
% For stats
cfg.keepindividual = ('yes');
cfg.method         = ('within');
% Baseline
grandavg_sta_baseline = ft_freqgrandaverage(cfg, StA_baseline{1}, StA_baseline{2}, StA_baseline{3}, StA_baseline{4}, StA_baseline{5}, StA_baseline{6}, StA_baseline{7}, StA_baseline{8}, StA_baseline{9}, StA_baseline{10}, StA_baseline{11}, StA_baseline{12}, StA_baseline{13}, StA_baseline{14}, StA_baseline{15}, StA_baseline{16}, StA_baseline{17}, StA_baseline{18}, StA_baseline{19}, StA_baseline{20},StA_baseline{21});
grandavg_cont_baseline = ft_freqgrandaverage(cfg, Cont_baseline{1}, Cont_baseline{2}, Cont_baseline{3}, Cont_baseline{4}, Cont_baseline{5}, Cont_baseline{6}, Cont_baseline{7}, Cont_baseline{8}, Cont_baseline{9}, Cont_baseline{10}, Cont_baseline{11}, Cont_baseline{12}, Cont_baseline{13}, Cont_baseline{14}, Cont_baseline{15}, Cont_baseline{16}, Cont_baseline{17}, Cont_baseline{18}, Cont_baseline{19}, Cont_baseline{20},Cont_baseline{21});
% Activation
grandavg_sta_activation = ft_freqgrandaverage(cfg, StA_activation{1}, StA_activation{2}, StA_activation{3}, StA_activation{4}, StA_activation{5}, StA_activation{6}, StA_activation{7}, StA_activation{8}, StA_activation{9}, StA_activation{10}, StA_activation{11}, StA_activation{12}, StA_activation{13}, StA_activation{14}, StA_activation{15}, StA_activation{16}, StA_activation{17}, StA_activation{18}, StA_activation{19}, StA_activation{20},StA_activation{21});
grandavg_cont_activation = ft_freqgrandaverage(cfg, Cont_activation{1}, Cont_activation{2}, Cont_activation{3}, Cont_activation{4}, Cont_activation{5}, Cont_activation{6}, Cont_activation{7}, Cont_activation{8}, Cont_activation{9}, Cont_activation{10}, Cont_activation{11}, Cont_activation{12}, Cont_activation{13}, Cont_activation{14}, Cont_activation{15}, Cont_activation{16}, Cont_activation{17}, Cont_activation{18}, Cont_activation{19}, Cont_activation{20},Cont_activation{21});

%% Time must be identical
% see http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#within-trial-experiments
    grandavg_sta_baseline.time = grandavg_sta_activation.time;
    grandavg_cont_baseline.time = grandavg_cont_activation.time;
%% REPMAT powspectrum
    grandavg_sta_baseline.powspctrm = repmat(grandavg_sta_baseline.powspctrm, [1,1,1,baseline_multiplier]);
    grandavg_cont_baseline.powspctrm = repmat(grandavg_cont_baseline.powspctrm, [1,1,1,baseline_multiplier]);
    
% Test is the time structure is correct: time_test = length(squeeze(grandavg_sta_baseline.powspctrm(1,1,1,:)))    
%% STATs 
% STA
cfg = [];
cfg.channel          = [1:64];
% cfg.latency          = [1 1.5]; % Just trying to test the alpha 100-500ms
cfg.method           = 'montecarlo';
cfg.frequency        = frequency;
cfg.statistic        = 'ft_statfun_actvsblT'; % 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
% Neighbours
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_activation{1});
cfg.avgoverfreq = 'yes';
% Groups
% Design
subj = 21; % put your number of subjects here
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

% Stats
[stats_StA_ep3_late_beta] = ft_freqstatistics(cfg, grandavg_sta_activation, grandavg_sta_baseline);
[stats_StA_ep2_early_alpha] = ft_freqstatistics(cfg, grandavg_sta_activation, grandavg_sta_baseline);

% save('stats_StA_ep3_late_beta','stats_StA_ep3_late_beta');
% save('stats_StA_ep2_early_alpha','stats_StA_ep2_early_alpha');

%                           STATs CONT
cfg = [];
cfg.channel          = [1:64];
% cfg.latency          = [1 1.5]; % Just trying to test the alpha 100-500ms
cfg.method           = 'montecarlo';
cfg.frequency        = frequency;
cfg.statistic        = 'ft_statfun_actvsblT'; % 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
% Neighbours
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, Cont_activation{1});
cfg.avgoverfreq = 'yes';
% Groups
% Design
subj = 21; % put your number of subjects here
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

% Stats
[stats_Cont_ep3_late_beta] = ft_freqstatistics(cfg, grandavg_cont_activation, grandavg_cont_baseline);
[stats_Cont_ep2_early_alpha] = ft_freqstatistics(cfg, grandavg_cont_activation, grandavg_cont_baseline);

% save('stats_Cont_ep3_late_beta','stats_Cont_ep3_late_beta');
% save('stats_Cont_ep2_early_alpha','stats_Cont_ep2_early_alpha');


%% STATS using single subject approach
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#within-trial-experiments

% Time must be identical
% see http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/#within-trial-experiments
for i = 1:21
    StA_baseline{i,1}.time = StA_activation{i,1}.time;
    Cont_baseline{i,1}.time = Cont_activation{i,1}.time;
end
% REPMAT powspectrum
for i = 1:21
    StA_baseline{i,1}.powspctrm = repmat(StA_baseline{i,1}.powspctrm, [1,1,baseline_multiplier]);
    Cont_baseline{i,1}.powspctrm = repmat(Cont_baseline{i,1}.powspctrm, [1,1,baseline_multiplier]);
end
for i =1
cfg = [];
cfg.channel          = [1:64];
% cfg.latency          = [0 0.1]; % Just trying to test the alpha 100-500ms
cfg.method           = 'montecarlo';
cfg.frequency        = frequency;
cfg.statistic        = 'ft_statfun_actvsblT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
% Neighbours
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_activation{i,1});
% Design
ntrials = size(StA_activation{i,1}.powspctrm,1);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];
cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
% Stats
[StA_stat] = ft_freqstatistics(cfg, StA_activation{i,1}, StA_baseline{i,1});
end
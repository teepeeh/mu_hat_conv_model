%% Plotting Conv TFR data

% PLots the time - frew imagsc plots for epsi2 and outcomes 
%% Paths
rmpath(genpath('E:\'));
clear all;
%% Directory
rootFolder = 'C:\Users\thoma\PhD\tosc\tph_convolution_modelling\tph_conv_mu_hat';
cd (rootFolder);
% Add Functions
addpath(genpath([rootFolder '\Functions\']));
addpath(genpath([rootFolder '\Data\']));
% Add SPM
spmpath = 'E:\spm12';
addpath(spmpath);
spm('defaults', 'eeg');
% Add EEGLAB
eeglabpath = 'E:\eeglab14_1_2b';
addpath(eeglabpath);

%% FIELD TRIP

% Loading in an old FT file as a template
% addpath E:\tpDATA\EEG_data\tanxiety_oscillation_analysis\Wavelet\tanx_osc_30_01_2020\Wavelet\Data;

% Toolboxes
ftpath = 'E:\fieldtrip-20181024';
% Add FieldTrip
addpath(genpath(ftpath));
ft_defaults;

load('tanx_wavelet_Cont_win_IQR_AB_baselined.mat')
ft_template = Cont_wavelet_win_AB{1, 1};

rmpath(genpath(ftpath));

%% Load each SPM file and save fttimelock to a FieldTrip template

% Groups
expgroup = [20 21 23:43];
expgroup = setdiff(expgroup, [36,37]);
contgroup = [1:19 22 37];
% Conditions
stimuli_blue_left   = 1;
stimuli_blue_right  = 2;
response_left       = 3;
response_right      = 4;
outcome_win         = 5;
outcome_lose        = 6;
regressor_muhat2    = 9; % Check


% State Anxiety
for i = expgroup
    % Load CONV MODEL continous TF data
    filetag = convertStringsToChars(num2str(i)); % nametag of SPM datafile
    subjfolder = convertStringsToChars(sprintf('participant_%d', i));
    D = spm_eeg_load([rootFolder '\Data\SPM_tf\Stimuli_locked' subjfolder '\Two_stimuli_regressor\conv_model_atf_d_' filetag '_tosc_spm_cont.mat']);
    % Stim left
    ft_sta_stim_left = ft_template;
    ft_sta_stim_left.time = D.fttimelock.time;
    ft_sta_stim_left.freq = D.fttimelock.freq;
    ft_sta_stim_left.powspctrm = [];
    ft_sta_stim_left.powspctrm = squeeze(D.fttimelock.powspctrm(stimuli_blue_left,:,:,:));
    ft_conv_sta_stim_left{i,1} = ft_sta_stim_left;
    clear ft_sta_stim_left
    % Stim right
    ft_sta_stim_right = ft_template;
    ft_sta_stim_right.time = D.fttimelock.time;
    ft_sta_stim_right.freq = D.fttimelock.freq;
    ft_sta_stim_right.powspctrm = [];
    ft_sta_stim_right.powspctrm = squeeze(D.fttimelock.powspctrm(stimuli_blue_right,:,:,:));
    ft_conv_sta_stim_right{i,1} = ft_sta_stim_right;
    clear ft_sta_stim_right
    % Resp left
    ft_sta_resp_left = ft_template;
    ft_sta_resp_left.time = D.fttimelock.time;
    ft_sta_resp_left.freq = D.fttimelock.freq;
    ft_sta_resp_left.powspctrm = [];
    ft_sta_resp_left.powspctrm = squeeze(D.fttimelock.powspctrm(response_left,:,:,:));
    ft_conv_sta_resp_left{i,1} = ft_sta_resp_left;
    clear ft_sta_resp_left
    % Resp right
    ft_sta_resp_right = ft_template;
    ft_sta_resp_right.time = D.fttimelock.time;
    ft_sta_resp_right.freq = D.fttimelock.freq;
    ft_sta_resp_right.powspctrm = [];
    ft_sta_resp_right.powspctrm = squeeze(D.fttimelock.powspctrm(response_right,:,:,:));
    ft_conv_sta_resp_right{i,1} = ft_sta_resp_right;
    clear ft_sta_resp_right    
    % Win
    ft_sta_win = ft_template;
    ft_sta_win.time = D.fttimelock.time;
    ft_sta_win.freq = D.fttimelock.freq;
    ft_sta_win.powspctrm = [];
    ft_sta_win.powspctrm = squeeze(D.fttimelock.powspctrm(outcome_win,:,:,:));
    ft_conv_sta_win{i,1} = ft_sta_win;
    clear ft_sta_win
    % Lose
    ft_sta_lose = ft_template;
    ft_sta_lose.time = D.fttimelock.time;
    ft_sta_lose.freq = D.fttimelock.freq;
    ft_sta_lose.powspctrm = [];
        ft_sta_lose.powspctrm = squeeze(D.fttimelock.powspctrm(outcome_lose,:,:,:));
        ft_conv_sta_lose{i,1} = ft_sta_lose;
    clear ft_sta_lose
    % MU HAT 2
    ft_sta_muhat2 = ft_template;
    ft_sta_muhat2.time = D.fttimelock.time;
    ft_sta_muhat2.freq = D.fttimelock.freq;
    ft_sta_muhat2.powspctrm = [];
    if length(D.conditions) == 7
        ft_sta_muhat2.powspctrm = squeeze(D.fttimelock.powspctrm(regressor_muhat2,:,:,:));
        ft_conv_sta_muhat2{i,1} = ft_sta_muhat2;
    else
        ft_sta_muhat2.powspctrm = squeeze(D.fttimelock.powspctrm(3,:,:,:));
        ft_conv_sta_muhat2{i,1} = ft_sta_muhat2;
    end
    clear ft_sta_epsi2
    clear S D
end
% For any missing participants
ft_conv_sta_stim_left = ft_conv_sta_stim_left(~cellfun('isempty',ft_conv_sta_stim_left));
ft_conv_sta_stim_right = ft_conv_sta_stim_right(~cellfun('isempty',ft_conv_sta_stim_right));
ft_conv_sta_resp_left = ft_conv_sta_resp_left(~cellfun('isempty',ft_conv_sta_resp_left));
ft_conv_sta_resp_right = ft_conv_sta_resp_right(~cellfun('isempty',ft_conv_sta_resp_right));
ft_conv_sta_win = ft_conv_sta_win(~cellfun('isempty',ft_conv_sta_win));
ft_conv_sta_lose = ft_conv_sta_lose(~cellfun('isempty',ft_conv_sta_lose));
ft_conv_sta_muhat2 = ft_conv_sta_muhat2(~cellfun('isempty',ft_conv_sta_muhat2));

% contgroup
for i = contgroup
    % Load CONV MODEL continous TF data
    filetag = convertStringsToChars(num2str(i)); % nametag of SPM datafile
    subjfolder = convertStringsToChars(sprintf('participant_%d', i));
    D = spm_eeg_load([rootFolder '\Data\SPM_tf\' subjfolder '\conv_model_atf_d_' filetag '_tosc_spm_cont.mat']);
    % Stim left
    ft_cont_stim_left = ft_template;
    ft_cont_stim_left.time = D.fttimelock.time;
    ft_cont_stim_left.freq = D.fttimelock.freq;
    ft_cont_stim_left.powspctrm = [];
    ft_cont_stim_left.powspctrm = squeeze(D.fttimelock.powspctrm(stimuli_blue_left,:,:,:));
    ft_conv_cont_stim_left{i,1} = ft_cont_stim_left;
    clear ft_cont_stim_left
    % Stim right
    ft_cont_stim_right = ft_template;
    ft_cont_stim_right.time = D.fttimelock.time;
    ft_cont_stim_right.freq = D.fttimelock.freq;
    ft_cont_stim_right.powspctrm = [];
    ft_cont_stim_right.powspctrm = squeeze(D.fttimelock.powspctrm(stimuli_blue_right,:,:,:));
    ft_conv_cont_stim_right{i,1} = ft_cont_stim_right;
    clear ft_cont_stim_right
    % Resp left
    ft_cont_resp_left = ft_template;
    ft_cont_resp_left.time = D.fttimelock.time;
    ft_cont_resp_left.freq = D.fttimelock.freq;
    ft_cont_resp_left.powspctrm = [];
    ft_cont_resp_left.powspctrm = squeeze(D.fttimelock.powspctrm(response_left,:,:,:));
    ft_conv_cont_resp_left{i,1} = ft_cont_resp_left;
    clear ft_cont_resp_left
    % Resp right
    ft_cont_resp_right = ft_template;
    ft_cont_resp_right.time = D.fttimelock.time;
    ft_cont_resp_right.freq = D.fttimelock.freq;
    ft_cont_resp_right.powspctrm = [];
    ft_cont_resp_right.powspctrm = squeeze(D.fttimelock.powspctrm(response_right,:,:,:));
    ft_conv_cont_resp_right{i,1} = ft_cont_resp_right;
    clear ft_cont_resp_right    
    % Win
    ft_cont_win = ft_template;
    ft_cont_win.time = D.fttimelock.time;
    ft_cont_win.freq = D.fttimelock.freq;
    ft_cont_win.powspctrm = [];
    ft_cont_win.powspctrm = squeeze(D.fttimelock.powspctrm(outcome_win,:,:,:));
    ft_conv_cont_win{i,1} = ft_cont_win;
    clear ft_cont_win
    % Lose
    ft_cont_lose = ft_template;
    ft_cont_lose.time = D.fttimelock.time;
    ft_cont_lose.freq = D.fttimelock.freq;
    ft_cont_lose.powspctrm = [];
        ft_cont_lose.powspctrm = squeeze(D.fttimelock.powspctrm(outcome_lose,:,:,:));
        ft_conv_cont_lose{i,1} = ft_cont_lose;
    clear ft_cont_lose
    % MU HAT 2
    ft_cont_muhat2 = ft_template;
    ft_cont_muhat2.time = D.fttimelock.time;
    ft_cont_muhat2.freq = D.fttimelock.freq;
    ft_cont_muhat2.powspctrm = [];
    if length(D.conditions) == 7
        ft_cont_muhat2.powspctrm = squeeze(D.fttimelock.powspctrm(regressor_muhat2,:,:,:));
        ft_conv_cont_muhat2{i,1} = ft_cont_muhat2;
    else
        ft_cont_muhat2.powspctrm = squeeze(D.fttimelock.powspctrm(3,:,:,:));
        ft_conv_cont_muhat2{i,1} = ft_cont_muhat2;
    end
    clear ft_cont_epsi2
    clear S D
end
% For any missing participants
ft_conv_cont_stim_left = ft_conv_cont_stim_left(~cellfun('isempty',ft_conv_cont_stim_left));
ft_conv_cont_stim_right = ft_conv_cont_stim_right(~cellfun('isempty',ft_conv_cont_stim_right));
ft_conv_cont_resp_left = ft_conv_cont_resp_left(~cellfun('isempty',ft_conv_cont_resp_left));
ft_conv_cont_resp_right = ft_conv_cont_resp_right(~cellfun('isempty',ft_conv_cont_resp_right));
ft_conv_cont_win = ft_conv_cont_win(~cellfun('isempty',ft_conv_cont_win));
ft_conv_cont_lose = ft_conv_cont_lose(~cellfun('isempty',ft_conv_cont_lose));
ft_conv_cont_muhat2 = ft_conv_cont_muhat2(~cellfun('isempty',ft_conv_cont_muhat2));

%% Load these fieldtrip convolution files into the grandavaerage

cfg = [];
cfg.channel = [1:64];
cfg.latency        = ('all');
% For plotting
cfg.keepindividual = ('no');
cfg.method         = ('within');

% Controls
grandavg_cont_stimBL = ft_freqgrandaverage(cfg, ft_conv_cont_stim_left{:});
grandavg_cont_stimBR = ft_freqgrandaverage(cfg, ft_conv_cont_stim_right{:});
grandavg_cont_respL = ft_freqgrandaverage(cfg, ft_conv_cont_resp_left{:});
grandavg_cont_respR = ft_freqgrandaverage(cfg, ft_conv_cont_resp_right{:});
grandavg_cont_win = ft_freqgrandaverage(cfg, ft_conv_cont_win{:});
grandavg_cont_lose = ft_freqgrandaverage(cfg, ft_conv_cont_lose{:});
grandavg_cont_muhat2 = ft_freqgrandaverage(cfg, ft_conv_cont_muhat2{:});
% State Anxiety
grandavg_sta_stimBL = ft_freqgrandaverage(cfg, ft_conv_sta_stim_left{:});
grandavg_sta_stimBR = ft_freqgrandaverage(cfg, ft_conv_sta_stim_right{:});
grandavg_sta_respL = ft_freqgrandaverage(cfg, ft_conv_sta_resp_left{:});
grandavg_sta_respR = ft_freqgrandaverage(cfg, ft_conv_sta_resp_right{:});
grandavg_sta_win = ft_freqgrandaverage(cfg, ft_conv_sta_win{:});
grandavg_sta_lose = ft_freqgrandaverage(cfg, ft_conv_sta_lose{:});
grandavg_sta_muhat2 = ft_freqgrandaverage(cfg, ft_conv_sta_muhat2{:});

%% Making DIFF Lose-win
contdiff = grandavg_cont_win;
contdiff.powspctrm = grandavg_cont_lose.powspctrm - grandavg_cont_win.powspctrm;

stadiff = grandavg_sta_win;
stadiff.powspctrm = grandavg_sta_lose.powspctrm - grandavg_sta_win.powspctrm;

%%                        plots within each group

ChanIndex = find(contains(grandavg_sta_stimBL.label,'FCz'));

%%                              Controls

close all;

% Stim Blue Left
subplot(5,2,1)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Stimulus Blue Left';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_stimBL.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_stimBL)

% Stim blue right
subplot(5,2,2)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Stimulus Blue Right';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_stimBR.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_stimBR)


% Resp Left
subplot(5,2,3)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Response Left';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_respL.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_respL)

% Resp Right
subplot(5,2,4)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Response Right';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_respR.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_respR)


% Win
subplot(5,2,5)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Win';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_win.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_win)

% Lose
subplot(5,2,6)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Lose';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_cont_lose.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_lose)


% Lose-Win
subplot(5,2,7:8)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Lose-Win';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-4, 4];
cfg.channel = contdiff.label(ChanIndex);
ft_singleplotTFR(cfg,contdiff)

% Epsi2
subplot(5,2,9:10)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'epsi2';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-200, 200];
cfg.channel = grandavg_cont_muhat2.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_cont_muhat2)

% Title
sgtitle('Controls')
% Background Colour
set(gcf,'color','white')


%%                          STATE ANXIETY
% StA
figure

% Stim Blue Left
subplot(5,2,1)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Stimulus Blue Left';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_stimBL.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_stimBL)

% Stim blue right
subplot(5,2,2)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Stimulus Blue Right';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_stimBR.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_stimBR)


% Resp Left
subplot(5,2,3)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Response Left';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_respL.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_respL)

% Resp Right
subplot(5,2,4)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Response Right';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_respR.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_respR)


% Win
subplot(5,2,5)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Win';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_win.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_win)

% Lose
subplot(5,2,6)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Lose';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-5, 5];
cfg.channel = grandavg_sta_lose.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_lose)


% Lose-Win
subplot(5,2,7:8)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'Lose-Win';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-4, 4];
cfg.channel = stadiff.label(ChanIndex);
ft_singleplotTFR(cfg,stadiff)

% Epsi2
subplot(5,2,9:10)
cfg = [];
% cfg.baseline = 'yes';
cfg.title = 'epsi2';
cfg.fontsize = 12;
cfg.zlim = 'maxabs'; % [-200, 200];
cfg.channel = grandavg_sta_muhat2.label(ChanIndex);
ft_singleplotTFR(cfg,grandavg_sta_muhat2)

% Title
sgtitle('State Anxiety')
% Background Colour
set(gcf,'color','white')

%%                              Between Groups

% epsi

% epsi2diff = grandavg_sta_muhat2;
% epsi2diff.powspctrm = grandavg_sta_muhat2.powspctrm - grandavg_cont_muhat2.powspctrm;
% 
% 
% % Plot
% 
% close all;
% cfg = [];
% % cfg.baseline = 'yes';
% cfg.title = 'epsi2';
% cfg.fontsize = 12;
% cfg.zlim = 'maxabs'; % [-5, 5];
% cfg.channel = epsi2diff.label(ChanIndex);
% ft_singleplotTFR(cfg,epsi2diff)
% % Title
% sgtitle('Between Group (StA-Cont)')
% % Background Colour
% set(gcf,'color','white')


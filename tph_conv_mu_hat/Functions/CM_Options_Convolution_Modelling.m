
%% Convolution Modelling - Settings Options

% Written by T. Hein June 2020
% Use this script to set all the options for convolution modelling
% functions included in 'tp_ConvModel.m'


% Updated for *** RESPONSE LOCKED 


%% Total Subj
options.subj_tot = [1];

%% Step1 - Events and Artefact Rejection

% Define events

options.events_artefacts.outcome_tp_events  = [140:142 240 241:242];       % Outcome  Locked = TB1 [140 win, 141, lose, 142, no resp] TB2 = [240, win, 241, lose, 242, no resp]
options.events_artefacts.stimuli_tp_events  = [110 111 210 211];           % Stimuli  Locked = TB1 [110 = blue left // 111 = blue right] TB2 [210 = blue left // 211 = blue right]
options.events_artefacts.response_tp_events = [121 124 127 221 224 227];   % Response Locked = TB1 = [121 124 127] left, right , no resp TB2 = [221 224 227] L, R, NR

% Artefact Rejection Indexes
options.events_artefacts.aRejindex   = ('GLM_aRej_osc_gamma_index_IQRx1_5.mat');
options.events_artefacts.aRej10index = ('EEGLAB_gamma_fix_10_aRej.mat');
options.events_artefacts.aRej14index = ('EEGLAB_gamma_fix_14_aRej.mat');
options.events_artefacts.aRej28index = ('EEGLAB_gamma_fix_28_aRej.mat');
options.events_artefacts.aRej30index = ('EEGLAB_gamma_fix_30_aRej.mat');
options.events_artefacts.aRej39index = ('EEGLAB_gamma_fix_39_aRej.mat');

%% Step2 - Convert 2 SPM and TF analysis
% Analysis Name
options.convertTF.SPM_prefix =  'SPM_tf\Stimuli_locked';
% Set the desired sampling freq for the TF amplitude data
options.sampling_freq = 256;
% Set time-freq analysis freq range
options.convertTF.tffreq = [8:2:30];
options.convertTF.channels = {'EEG'};
options.convertTF.timewin = [-Inf Inf];
options.convertTF.method = 'morlet';
options.convertTF.settings.ncycles = 5;

... cfg.width      = 7;  % GAMMA [30:2:100]
... cfg.width      = 5;  % BETA  [8:2:30]
... cfg.width      = 3;  % ALPHA [8:2:30]
... cfg.width      = 2;  % THETA [4:1:8]

%% Convolution Modelling
% Data DIR
options.conv.data_dir = [rootFolder '\Data\SPM_tf\'];
% Convolution Modelling time window of analysis
options.conv.toi = [-200 2000];
% Sensor analysis (EEG for all sensors, or e.g. Cz for one)
options.conv.sensors = 'EEG';
% Order of Fourier-Basis (how many sine & cosine functions)
options.conv.ofb = 12;

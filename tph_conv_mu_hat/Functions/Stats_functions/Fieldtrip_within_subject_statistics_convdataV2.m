%% #2 Within group tests FIELDTRIP on convolution data

% Notes:

% Below Maria updated the code for both finding the TOI using the sampling
% frequency and using milliseconds and an index. Using the sampling
% frequency is now obsoltete, but below is the updated code for future
% reference.

% Find TOI using sampling freq
% sampling_freq = 256; % sampling points
% baseline_duration = 200; % in ms
% baseline_nsampl = ceil((sampling_freq*baseline_duration)/1000);

%% Remove paths
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
% FIELD TRIP
ftpath = 'E:\fieldtrip-20181024';
%% Make a Fieldtrip template for SPM data
load('tanx_wavelet_Cont_win_IQR_AB_baselined.mat')
ft_template = Cont_wavelet_win_AB{1, 1};
%% Load each SPM file and save fttimelock to a FieldTrip template
% Exp groups
expgroup = [20 21 23:43];
expgroup = setdiff(expgroup, [36,37]);
contgroup = [1:19 22 37];
% Convolution conditions
condition_win = 1;
condition_lose = 2;
condition_noresp = 3;
condition_epsi2 = 4;
condition_epsi3 = 5;
%% Choose regressor + time window
% Epsilon 2 cases
case1 = 'epsilon2 early alpha';
case4 = 'epsilon2 late alpha';
case5 = 'epsilon2 late beta';
% Epsilon 3 cases
case2 = 'epsilon3 late beta';
case3 = 'epsilon3 early alpha';
case6 = 'epsilon3 late alpha';
% For loop
cases = {'epsilon2 early alpha';'epsilon3 late beta';'epsilon3 early alpha';'epsilon2 late alpha';'epsilon2 late beta';'epsilon3 late alpha'};
%% RUN
for k = 1:6
% mycase = case1;
mycase = cell2mat(cases(k));
% Remove FT
rmpath(genpath(ftpath));
% add SPM 
addpath(spmpath);
spm('defaults', 'eeg');
% Run
switch mycase
    case 'epsilon2 early alpha'
        % STATE ANXIETY group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_sta_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_sta_epsi2 = ft_conv_sta_epsi2(~cellfun('isempty',ft_conv_sta_epsi2));
        % CONTROL group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_cont_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_cont_epsi2 = ft_conv_cont_epsi2(~cellfun('isempty',ft_conv_cont_epsi2));
    case 'epsilon2 late alpha'
        % STATE ANXIETY group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_sta_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_sta_epsi2 = ft_conv_sta_epsi2(~cellfun('isempty',ft_conv_sta_epsi2));
        % CONTROL group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_cont_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_cont_epsi2 = ft_conv_cont_epsi2(~cellfun('isempty',ft_conv_cont_epsi2));
    case 'epsilon2 late beta'
        % STATE ANXIETY group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_sta_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_sta_epsi2 = ft_conv_sta_epsi2(~cellfun('isempty',ft_conv_sta_epsi2));
        % CONTROL group
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
            ft_template.powspctrm = squeeze(D.fttimelock.powspctrm(condition_epsi2,:,:,:));
            ft_conv_cont_epsi2{i,1} = ft_template;
            clear S D
        end
        % For any missing participants
        ft_conv_cont_epsi2 = ft_conv_cont_epsi2(~cellfun('isempty',ft_conv_cont_epsi2));
    case 'epsilon3 late beta'
        % STATE ANXIETY group
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
        % CONTROL group
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
    case 'epsilon3 early alpha'
        % STATE ANXIETY group
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
        % CONTROL group
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
    case 'epsilon3 late alpha'
        % STATE ANXIETY group
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
        % CONTROL group
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
end
%% Update FT Path
rmpath(spmpath);
addpath(genpath(ftpath));
ft_defaults;
%% Extracting data based on TOI = new powspctrm + time

switch mycase
    case 'epsilon2 early alpha'
        frequency = [8 12];
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi2{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi2{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi2{i};
            Cont_bl{i,1}.time = Cont_bl{i,1}.time([bl_start:bl_end]);
            Cont_bl{i,1}.powspctrm = Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [0.1 0.5]; %in seconds, careful, fieldtrip time field is in seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi2{i};
            StA_act{i,1}.time = ft_conv_sta_epsi2{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi2{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi2{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
        
    case 'epsilon2 late alpha'
        frequency = [8 12];
        bl_repmat = 3;
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi2{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi2{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi2{i};
            Cont_bl{i,1}.time = Cont_bl{i,1}.time([bl_start:bl_end]);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [1.1 1.6]; %in seconds, careful, fieldtrip time field is in seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi2{i};
            StA_act{i,1}.time = ft_conv_sta_epsi2{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi2{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi2{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
        
    case 'epsilon2 late beta'
        frequency = [13 25];
        bl_repmat = 3;
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi2{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi2{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi2{i};
            Cont_bl{i,1}.time = Cont_bl{i,1}.time([bl_start:bl_end]);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [1.0 1.6]; %in seconds, careful, fieldtrip time field is in seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi2{i};
            StA_act{i,1}.time = ft_conv_sta_epsi2{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi2{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi2{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
        
    case 'epsilon3 late beta'
        frequency = [13 25];
        bl_repmat = 2;
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi3{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi3{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi3{i};
            Cont_bl{i,1}.time = Cont_bl{i,1}.time([bl_start:bl_end]);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [1.4 1.7]; % Seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi3{i};
            StA_act{i,1}.time = ft_conv_sta_epsi3{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi3{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi3{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
        
        
    case 'epsilon3 early alpha'
        frequency = [8 12];
        bl_repmat = 2;
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi3{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi3{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi3{i};
            Cont_bl{i,1}.time = Cont_bl{i,1}.time([bl_start:bl_end]);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [0.1 0.5]; % Seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi3{i};
            StA_act{i,1}.time = ft_conv_sta_epsi3{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi3{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi3{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
        
    case 'epsilon3 late alpha'
        frequency = [8 12];
        bl_repmat = 0;
        % Baseline
        % TOI using time
        mytime = ft_conv_sta_epsi3{1}.time; % time vector in ms  -200:2000 ms in one participant
        blTOI = [-0.2 0]; % Seconds
        bl_start = find(min(abs(blTOI(1) - mytime)) == abs(blTOI(1)  - mytime));
        bl_end = find(min(abs(blTOI(2) - mytime)) == abs(blTOI(2)  - mytime));
        % Extract data
        for i = 1:21
            % State anxiety
            StA_bl{i,1} = ft_conv_sta_epsi3{i};
            StA_bl{i,1}.time  = StA_bl{i,1}.time([bl_start:bl_end]);
            StA_bl{i,1}.powspctrm = nanmean(StA_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
            % Controls
            Cont_bl{i,1} = ft_conv_cont_epsi3{i};
            Cont_bl{i,1}.time = nanmean(Cont_bl{i,1}.time([bl_start:bl_end]),3);
            Cont_bl{i,1}.powspctrm = nanmean(Cont_bl{i,1}.powspctrm(:,:,bl_start:bl_end),3);
        end
        
        % Activation
        actTOI = [1.5 1.7]; % Seconds
        act_start = find(min(abs(actTOI(1) - mytime)) == abs(actTOI(1)  - mytime));
        act_end = find(min(abs(actTOI(2) - mytime)) == abs(actTOI(2)  - mytime));
        for i = 1:21
            % State anxiety
            StA_act{i,1} = ft_conv_sta_epsi3{i};
            StA_act{i,1}.time = ft_conv_sta_epsi3{i,1}.time([act_start:act_end]);
            StA_act{i,1}.powspctrm = StA_act{i,1}.powspctrm(:,:,[act_start:act_end]);
            % Controls
            Cont_act{i,1} = ft_conv_cont_epsi3{i};
            Cont_act{i,1}.time = ft_conv_cont_epsi3{i,1}.time([act_start:act_end]);
            Cont_act{i,1}.powspctrm = Cont_act{i,1}.powspctrm(:,:,[act_start:act_end]);
        end
        % Repmat Baseline
        for i = 1:21
            base_repmat_sta = (StA_bl{i,1}.powspctrm);
            base_repmat_cont = (Cont_bl{i,1}.powspctrm);
            for n = 2:length([act_start:act_end])
                StA_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_sta, [1,1]);
                Cont_bl{i,1}.powspctrm(:,:,n) = repmat(base_repmat_cont, [1,1]);
            end
            clear base_repmat_sta base_repmat_cont
        end
        
        % Baseline time
        for i = 1:21
            StA_bl{i,1}.time = StA_act{i,1}.time;
            Cont_bl{i,1}.time = Cont_act{i,1}.time;
        end
end

%% FT grand avaerage
% Keeping individuals
cfg = [];
cfg.channel = [1:64];
cfg.latency        = ('all');
cfg.keepindividual = ('yes');
cfg.method         = ('within');
% Baseline
grandavg_sta_bl = ft_freqgrandaverage(cfg, StA_bl{:});
grandavg_cont_bl = ft_freqgrandaverage(cfg, Cont_bl{:});
% Activation
grandavg_sta_act = ft_freqgrandaverage(cfg, StA_act{:});
grandavg_cont_act = ft_freqgrandaverage(cfg, Cont_act{:});
%% STATS

stats_svdir = [rootFolder + "\Stats\within_baseline_fieldtrip"];

switch mycase
    case 'epsilon2 early alpha'
        cd ([stats_svdir + "\epsi2"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep2_early_alpha] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep2_early_alpha','stats_StA_ep2_early_alpha');   
        % CONT
        [stats_Cont_ep2_early_alpha] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep2_early_alpha','stats_Cont_ep2_early_alpha');
        
    case 'epsilon2 late alpha'
        cd ([stats_svdir + "\epsi2"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep2_late_alpha] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep2_late_alpha','stats_StA_ep2_late_alpha');
        % CONT
        [stats_Cont_ep2_late_alpha] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep2_late_alpha','stats_Cont_ep2_late_alpha');
        
    case 'epsilon2 late beta'
        cd ([stats_svdir + "\epsi2"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep2_late_beta] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep2_late_beta','stats_StA_ep2_late_beta');
        % CONT
        [stats_Cont_ep2_late_beta] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep2_late_beta','stats_Cont_ep2_late_beta');
        
    case 'epsilon3 late beta'
        cd ([stats_svdir + "\epsi3"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep3_late_beta] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep3_late_beta','stats_StA_ep3_late_beta');
        % CONT
        [stats_Cont_ep3_late_beta] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep3_late_beta','stats_Cont_ep3_late_beta');
        
    case 'epsilon3 early alpha'
        cd ([stats_svdir + "\epsi3"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep3_early_alpha] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep3_early_alpha','stats_StA_ep3_early_alpha'); 
        % CONT
        [stats_Cont_ep3_early_alpha] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep3_early_alpha','stats_Cont_ep3_early_alpha');
        
    case 'epsilon3 late alpha'
        cd ([stats_svdir + "\epsi3"]);
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
        cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, StA_act{1});
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
        [stats_StA_ep3_late_alpha] = ft_freqstatistics(cfg, grandavg_sta_act, grandavg_sta_bl);
%         save('stats_StA_ep3_late_alpha','stats_StA_ep3_late_alpha');
        % CONT
        [stats_Cont_ep3_late_alpha] = ft_freqstatistics(cfg, grandavg_cont_act, grandavg_cont_bl);
%         save('stats_Cont_ep3_late_alpha','stats_Cont_ep3_late_alpha');
end

clear grandavg_sta_act grandavg_sta_bl grandavg_cont_bl grandavg_cont_act...
    StA_bl Cont_bl StA_act Cont_act frequency mytime blTOI bl_start bl_end ...
    actTOI act_start act_end StA_act Cont_act design filetag
                  
end





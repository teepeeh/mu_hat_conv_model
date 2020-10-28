function IRs = CM_step3_Impulse_Response_Matrix(rootFolder, options)

%% Convolution Modelling Time-Frequency - Step 3: creating impulase response (wrapper)
% Written by T. Hein June 2020
% Load EEG files and artefact rejection and downsamples to create a matric
% of impulase responses and computational model regressor responses from
% the 'SPM_tosc_conttinous' files and behavioural data.


%% IR matrix

% Design matrix (pre-convolution), in which each column codes one outcome type / condition of interest:
% Load artefact rejection index
load(options.events_artefacts.aRejindex);
% Starts EEGLAB
eeglab; close all;

%% SAVE

% DIR for impulse response matrices
sample_rate_dir = (['sr', num2str(options.sampling_freq)]);
svdir = char([rootFolder + "\Data\IR_designmatrix"]);
svdir = fullfile(svdir, sample_rate_dir);
% Creates a folder if necessary
if ~exist([svdir],'dir')
    mkdir([svdir]);
end

cd(svdir);

%% RUN

for n = options.subj_tot
    
    %%                              Stimuli
    
    eeg_data_load = 'C:\Users\thoma\PhD\tosc\tph_convolution_modelling\tph_conv_mu_hat\Data\SPM_stimuli_locked_cont';
    addpath(eeg_data_load)
    % Input file name
    fn = sprintf('%d_tosc_stimuli_locked_aRej.set',n);
    % Load EEG file
    EEG = pop_loadset(fn);
    % Downsample to same as SPM final conv modelling
    EEG = pop_resample(EEG, options.sampling_freq);
    % Create matrix for impulse responses
    
    % USING ONE REGRESSOR
    %     stimuli = zeros(size(EEG.data,2),1,'single'); % COULD BE 2
    %     stimuli_key = {'Stimu Blue_R&L'}; % {'blue left';'blue right'};
    
    % USING TWO REGRESSORS
    stimuli = zeros(size(EEG.data,2),2,'single');
    stimuli_key = {'Stim_blue_left','Stim_blue_right'}; % {'blue left';'blue right'};
    
    % Create Matrix for HGF regressors
    HGF =  zeros(size(EEG.data,2),1,'single');
    HGF_key = {'|muhat2|'};
    % Getting event lables
    tpevents =  {EEG.event.type}';
    tpevents = cellfun(@str2double, tpevents);
    % Getting event latencies
    tplats =  {EEG.event.latency}';
    tplats = cell2mat(tplats);
    
    % Getting regressors subject's 3-lev HGF
    % Load behavioural data
    fname = sprintf('tanx_%d.mat',n);
    load (fname);
    r = player_struct.tapas.est;
    clear MUhat2;
    % HGF - epsi2
    MUhat2(:,1) = abs(r.traj.muhat(:,2)); %absolute
    % Here need to add in deleting the aRej
    aRej = GLM_aRej_osc_gamma_index_IQRx1_5{1, n};
    MUhat2(aRej) = [];
    
    % RUN outcomes/regressors according to latencies
    
    % Locked TB1 [110 = blue left // 111 = blue right] TB2 [210 = blue left // 211 = blue right]
    
    % For each trial take the event and time
    % USING ONE REGRESSOR
    
    %     for i = 1:size(tpevents,1)
    %         if isequal(tpevents(i),110) || isequal(tpevents(i),111)
    %             lat = tplats(i);
    %             stimuli(round(lat),1) = 1;
    %             HGF(round(lat),1) = MUhat2(i);
    %         elseif isequal(tpevents(i),210) || isequal(tpevents(i),211)
    %             lat = tplats(i);
    %             stimuli(round(lat),1) = 1;
    %             HGF(round(lat),1) = MUhat2(i);
    %         else
    %             lat = tplats(i);
    %             stimuli(round(lat),1) = 1;
    %             HGF(round(lat),1) = MUhat2(i);
    %         end
    %     end
    
    % USING TWO REGRESSORS
    % For each trial take the event and time
    for i = 1:size(tpevents,1)
        if isequal(tpevents(i),110) || isequal(tpevents(i),210)
            lat = tplats(i);
            stimuli(round(lat),1) = 1;
            HGF(round(lat),1) = MUhat2(i);
        elseif isequal(tpevents(i),111) || isequal(tpevents(i),211)
            lat = tplats(i);
            stimuli(round(lat),2) = 1;
            HGF(round(lat),1) = MUhat2(i);
        end
    end
    
    % Compiling
    IR_wrapper.IR_designmatrix.stimuli.stimuli{n} = stimuli;
    IR_wrapper.IR_designmatrix.stimuli.key = stimuli_key;
    IR_wrapper.HGF_designmatrix.stimuli.HGF{n} = HGF;
    IR_wrapper.HGF_designmatrix.stimuli.key = HGF_key;
    
    clear fn EEG lat HGF tpevents tplats
    %%                       Outcome locked
    
    rmpath(eeg_data_load);
    eeg_data_load = 'C:\Users\thoma\PhD\tosc\tph_convolution_modelling\tph_conv_mu_hat\Data\SPM_outcome_locked_cont';
    addpath(eeg_data_load)
    % Input file name
    fn = sprintf('%d_tosc_outcome_locked_aRej.set',n);
    % Load EEG file
    EEG = pop_loadset(fn);
    % Downsample to same as SPM final conv modelling
    EEG = pop_resample(EEG, options.sampling_freq);
    % Create matrix for impulse responses
    outcomes = zeros(size(EEG.data,2),3,'single');
    outcomes_key = {'win';'lose';'outcome: no response'};
    % Getting event lables
    tpevents =  {EEG.event.type}';
    tpevents = cellfun(@str2double, tpevents);
    % Getting event latencies
    tplats =  {EEG.event.latency}';
    tplats = cell2mat(tplats);
    
    % For each trial take the event and time
    for i = 1:size(tpevents,1)
        if isequal(tpevents(i),140) || isequal(tpevents(i),240)
            lat = tplats(i);
            outcomes(round(lat),1) = 1;
        elseif isequal(tpevents(i),141) || isequal(tpevents(i),241)
            lat = tplats(i);
            outcomes(round(lat),2) = 1;
        else
            lat = tplats(i);
            outcomes(round(lat),3) = 1;
        end
    end
    % Compiling
    IR_wrapper.IR_designmatrix.outcomes.outcomes{n} = outcomes;
    IR_wrapper.IR_designmatrix.outcomes.key = outcomes_key;
    
    clear fn EEG lat HGF tpevents tplats
    %%                        RESPONSE Locked
    rmpath(eeg_data_load);
    eeg_data_load = 'C:\Users\thoma\PhD\tosc\tph_convolution_modelling\tph_conv_mu_hat\Data\SPM_response_locked_cont';
    addpath(eeg_data_load)
    % Input file name
    fn = sprintf('%d_tosc_response_locked_aRej.set',n);
    % Load EEG file
    EEG = pop_loadset(fn);
    % Downsample to same as SPM final conv modelling
    EEG = pop_resample(EEG, options.sampling_freq);
    % Create matrix for impulse responses
    response = zeros(size(EEG.data,2),3,'single');
    response_key = {'left resp';'right resp';'response: no resp'};
    % Getting event lables
    tpevents =  {EEG.event.type}';
    tpevents = cellfun(@str2double, tpevents);
    % Getting event latencies
    tplats =  {EEG.event.latency}';
    tplats = cell2mat(tplats);
    
    % For each trial take the event and time
    for i = 1:size(tpevents,1)
        if isequal(tpevents(i),121) || isequal(tpevents(i),221) % Left
            lat = tplats(i);
            response(round(lat),1) = 1;
        elseif isequal(tpevents(i),124) || isequal(tpevents(i),224) % Right
            lat = tplats(i);
            response(round(lat),2) = 1;
        else
            lat = tplats(i);
            response(round(lat),3) = 1;
        end
    end
    % Compiling
    IR_wrapper.IR_designmatrix.response.response{n} = response;
    IR_wrapper.IR_designmatrix.response.key = response_key;
    %% Saving
    save([svdir + "\IR_wrapper_" + sample_rate_dir + ".mat"] ,'IR_wrapper');
    % Clearing
    clear EEG tpevents lat outcomes tplats fn HGF response stimuli
end

end

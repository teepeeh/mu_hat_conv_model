%% Maria's code for checking event values in EEG data (IR matrix issue)

%%
rmpath(genpath('C:\Users\thoma\PhD\tosc\tph_convolution_modelling\'));
clear all;
%% Directory
rootFolder = 'C:\Users\thoma\PhD\tosc\tph_convolution_modelling\tph_conv_mu_hat';
cd (rootFolder);
% Add Functions
addpath(genpath([rootFolder '\Functions\']));
addpath(genpath([rootFolder '\Data\']));
addpath(genpath([rootFolder '\Plots\']));
addpath(genpath([rootFolder '\Stats\']));
% Add EEGLAB
eeglabpath = 'C:\Users\thoma\PhD\eeglab14_1_2b';
addpath(eeglabpath);

%% Maria Code

EEGdata = '1_I64efMerged.set';
eeglab; close all;
EEG = pop_loadset(EEGdata);

srate = EEG.srate;
% evtpye = cell2struct({EEG.event.type}) % This was Maria's but didnt work
% for me??
evtpye = {EEG.event.type};
evlat = [EEG.event.latency];
evtime = evlat*1000/srate; %milliseconds
lengthcheck = (max(evtime)/1000)/60; % 26 minutes


options.events_artefacts.outcome_tp_events = [140:142 240:242];
% Outcome ?Locked

options.events_artefacts.stimuli_tp_events ? = [110 111 210 211];
% Stimuli ?Locked TB1 [110 = blue left // 111 = blue right] TB2 [210 = blue left // 211 = blue right]

options.events_artefacts.response_tp_events = [121 124 127 221 224 227];
% Response Locked TB1 = [121 124 127] left, right , no resp TB2 = [221 224 227] L, R, NR


indS= sort(cat(2, find(strcmpi(evtype, '110')), find(strcmpi(evtype, '111')), find(strcmpi(evtype, '210')),...
    find(strcmpi(evtype, '211'))));


indR= sort(cat(2, find(strcmpi(evtype, '121')), find(strcmpi(evtype, '124')), find(strcmpi(evtype, '127')),...
    find(strcmpi(evtype, '221')), find(strcmpi(evtype, '224')), find(strcmpi(evtype, '227'))));


indO= sort(cat(2, find(strcmpi(evtype, '140')), find(strcmpi(evtype, '141')), find(strcmpi(evtype, '142')),...
    find(strcmpi(evtype, '240')), find(strcmpi(evtype, '241')), find(strcmpi(evtype, '242'))));


%% time differences between Stimulus and Response: Reaction Time, should be
%%in range 200 - 600 ms
timediff_SR = evtime(indR) - evtime(indS);

figure;hist(timediff_SR,100)  %OK



%% time differences between Response and Outcome:
%%in 1000 - 1045 ms
timediff_RO = evtime(indO) - evtime(indR);
figure;hist(timediff_RO,100)  %OK










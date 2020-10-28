%% Convolution Modelling - Time-Frequency EEG data.
% Written by T. Hein June 2020
% Requires preprocessed EEG data + index of trials to be rejected for each
% subject.

% RESPONSE LOCKED 

%% Run 
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
% Add SPM
spmpath = 'C:\Users\thoma\PhD\spm12';
addpath(spmpath);
% Add EEGLAB
eeglabpath = 'C:\Users\thoma\PhD\eeglab14_1_2b';
addpath(eeglabpath);
%% Options
% Edit this script to set options for the following functions
CM_Options_Convolution_Modelling;
%% Step 1 - EGGLAB event selection + artefact rejection 
% Calls function
eEEGaRej = CM_step1_EEG_events_artefacts(rootFolder, options);
clearvars -except eEEGaRej rootFolder spmpath  eeglabpath options
clc;
% Check saved EEG files
dir ([rootFolder '\Data\SPM_tosc_continous\']);
%% Step 2 - Convert to SPM + time frequency analysis 
% If eEEGaRej (names of files for step2 function generated in step 1) is missing.
if exist('eEEGaRej','var') == 0    
    eEEGaRej_script;
end
% Calls function
spmTF = CM_step2_EEG2SPM_TF(rootFolder, eEEGaRej, options);
clearvars -except rootFolder spmpath spmTF eEEGaRej eeglabpath aRejindex options
% Check saved EEG files
dir ([rootFolder '\Data\SPM_tf\']);
cd(rootFolder);
%% Step 3 - Impulse Response Matrix
CM_step3_Impulse_Response_Matrix(rootFolder, options);
clearvars -except rootFolder spmpath spmTF eEEGaRej eeglabpath aRejindex options
cd(rootFolder);
clc;
%% Step 4 - Convolution Modelling
CM_step4_Convolution_Modelling(rootFolder, options);



%%                          UNUSED STEPS

%% Step 5 - Convert to SPM images + Smoothing
% CM_step5_SPM_images(rootFolder, options);
%% Stats - 2nd level one sample ttest (within)
% This will only work if all given subjects in options for stats have been 
% processed from step1-5.
% CM_2ndLevel_groupmean(rootFolder, options);
% close all;
%% Stats - 2nd level two sample ttest (between groups)
% This will only work if all given subjects in options for stats have been 
% processed from step1-5.
% CM_2ndLevel_groupDiff(rootFolder, options);
% close all;


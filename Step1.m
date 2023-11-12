% Thesis Sanna Bakels, 4480279
% Step 1: Filtering the raw EEG data and creating SSR

clear; close all; clc;

% Using the Matlab Toolbox EEGLAB from:
% Delorme, Arnaud, and Scott Makeig. "EEGLAB: an open source toolbox for 
% analysis of single-trial EEG dynamics including independent component 
% analysis." Journal of neuroscience methods 134.1 (2004): 9-21. 

%% Open the EEG file
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

EEG = pop_loadset('filename', 'EEG_subject2', 'filepath', '~/Matlab/ExperimentData/Subject 2/');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
EEG = eeg_checkset( EEG );
eeglab redraw;

%% Extract event 2 data 
EEG = pop_rmdat( EEG, {'  2'},[0 12.5] ,0);
EEG.setname='EEG trigger 2 events';
eeglab redraw;

%% Baseline removal + filtering
EEG = pop_rmbase( EEG, [],[]);
EEG.setname='baseline removed';

% Filter the data between 0.8 Hz and 120 Hz as well as removing the line noise
EEG = pop_eegfiltnew(EEG, 'locutoff',0.8,'hicutoff',120);
EEG.setname='filter 1';
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:64] ,'computepower',1,'linefreqs',[50 100] ,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
EEG = pop_reref( EEG, []);
EEG.setname='re-referenced';
eeglab redraw;

%% Manually reject bad channels (most of the times channel 13 and 19)
figure; 
pop_spectopo(EEG, 1, [0      250009.2773], 'EEG' , 'freq', [6 10 22], 'freqrange',[2 25],'electrodes','off');
eeglab redraw;

%% Apply ICA decomposition to clean the data
pop_processMARA ( ALLEEG,EEG,CURRENTSET )
eeglab redraw;

EEG.setname='ICA';
EEG = pop_iclabel(EEG, 'default');
EEG = pop_subcomp( EEG, [], 0);
EEG = eeg_checkset( EEG );
eeglab redraw;

%% Rereferce to average
EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, []);
eeglab redraw;

%% Compute the SSR for one subject, each trial, such that you have 160 periods in total
% Based on equation 1. from Vlaar et al. (2017) - Quantification of
% task-dependent cortical activation evoked by robotic continous wrist
% joint manipulation in chronic hemiparetic stroke. 

time = EEG.times;
numTimePoints = numel(time);
numTrials = 20;
numEpochs = EEG.trials;
numPeriods = numTrials * numEpochs;
numChannels = EEG.nbchan;

allData = zeros(numChannels, numTimePoints, numPeriods);
averageData = zeros(numChannels, numTimePoints, 1);
baseFileName = 'EEG_trial';

% Loop through and load each file one by one
for i = 1:numTrials
    filePath = fullfile(inputFolderPath, sprintf('%s%d.mat', baseFileName, i));
    load(sprintf(filePath));

    for j = 1:numEpochs
        % Calculate the correct index to store the data
        index = (i-1) * numEpochs + j;
        allData(:, :, index) = EEG.data(:,:,j);
    end
end

% Calculate the mean SSR for all channels and time points
averageData = mean(allData, 3); % Compute the average along the third dimension (epochs)
reshapedData = reshape(allData(1, :, :), [numTimePoints, numPeriods]);

subplot(1, 2, 1);
plot(time, reshapedData(:,1), 'LineWidth', 0.005, 'DisplayName', 'Electrical activity for one epoch');
hold on
plot(time, reshapedData(:,2:end), 'LineWidth', 0.005, 'HandleVisibility', 'off');
plot(time, averageData(1, :), 'k', 'LineWidth', 3, 'DisplayName', 'SSR')

subplot(1, 2, 2);
plot(time, averageData(1, :), 'k', 'LineWidth', 2);
function decoderFileLocationAndNames = calibrateMainAndSecondaryKalmanVelocityDecodersForBci(nevFilebase, nevFilesForTrain, trainParams, subject)
% Function operates under assumption that two decoders will be trained
% during calibration

global params codes

addpath(genpath('C:\Users\rigmdata\nevutils\'))
addpath(genpath('C:\Users\rigmdata\spikesort\'))
addpath(genpath('C:\Users\rigmdata\bciCode\'))
addpath(genpath('C:\Users\rigmdata\Documents\Ex\'))

netFolder = params.nasNetFolder;
nasNetName = trainParams.nasNetwork;
% Retrieve decoder parameters for both decoders
firstGamma = trainParams.firstGamma;
firstChannelNumbersUse = trainParams.firstRippleChannelNumbersInBci;
firstQvalue = trainParams.firstKalmanQ;
firstNumLatents = trainParams.firstNumberFaLatents;
firstVelocityToCalWith = trainParams.firstVelocityToCalibrateWith;

secondGamma = trainParams.secondGamma;
secondChannelNumbersUse = trainParams.secondRippleChannelNumbersInBci;
secondQvalue = trainParams.secondKalmanQ;
secondNumLatents = trainParams.secondNumberFaLatents;
secondVelocityToCalWith = trainParams.secondVelocityToCalibrateWith;

%% Reading in the NEV File
% don't double read the base if it's also a train file
if any(strcmp(nevFilesForTrain, nevFilebase))
    includeBaseForTrain = true;
else
    includeBaseForTrain = false;
end
nevFilesForTrain(strcmp(nevFilesForTrain, nevFilebase)) = [];
% Read in nevFileBase as a nev file and stores waveforms and base that
% contains channels, codes, and timestamps. NevFileBase is needed for the
% global parameters.
[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
% Reads in NEVS that contains training data that will be used for training
% decoder. 
for nevFlInd = 1:length(nevFilesForTrain)
    [nevNx, wavesNx] = readNEV(nevFilesForTrain{nevFlInd});
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    % Creates nev and waves variables that are concatenated across all NEVs
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end

% Convert global nev file to dat for use with matlab
datBase = nev2dat(nevBase, 'nevreadflag', 1);

% Retrieve network's labels that have not yet been thresholded yet.
% NevLabelledData is similar to nev and has 3 columns
[~,~, netLabels] = runNASNet({nev, waves},firstGamma, ...
    'netFolder', netFolder, 'netname', nasNetName);
% Generates Save Folder directory for decoders
bciDecoderSaveDrive = 'X:\';
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end

% For decoder filenames
timestamp = datestr(now, 'HH-MM-SS');
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));

% Train decoders 
fprintf('\n Training first decoder \n')
[M0, M1, M2, channelsKeep, A, Q, C, R, beta, K] = trainKalmanDecoder(nev, trainParams, params, codes, datBase, nevBase, ...
    netLabels, firstGamma, firstChannelNumbersUse, firstQvalue, firstNumLatents, firstVelocityToCalWith, includeBaseForTrain);

bciDecoderSaveName = sprintf('%s%sKalmanBciDecoderOne_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), timestamp);
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2', 'channelsKeep', 'A', 'Q', 'C', 'R', 'beta', 'K', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
firstDecoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);

fprintf('\n Training second decoder \n')
[M0, M1, M2, channelsKeep, A, Q, C, R, beta, K] = trainKalmanDecoder(nev, trainParams, params, codes, datBase, nevBase, ...
    netLabels, secondGamma, secondChannelNumbersUse, secondQvalue, secondNumLatents, secondVelocityToCalWith, includeBaseForTrain);

bciDecoderSaveName = sprintf('%s%sKalmanBciDecoderTwo_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), timestamp);
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2', 'channelsKeep', 'A', 'Q', 'C', 'R', 'beta', 'K', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
secondDecoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);

decoderFileLocationAndNames = strjoin({firstDecoderFileLocationAndName, secondDecoderFileLocationAndName}, '\n');

end
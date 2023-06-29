function [decoderFileLocationAndName] = trainDecoderFromNev(nevFilebase, subject)
% function [M0, M1, M2, channelsKeep, A, Q, C, R, beta, K] = trainKalmanDecoderFromNev(nev, trainParams, params, codes, datBase, nevBase,...
%         netLabels, gamma, channelNumbersUse, Qvalue, numLatents, velocityToCalibrateWith, includeBaseForTrain)

exGlobals;
dataPath = 'E:\';
xmlFile = 'bci_rewardAxisDecoder.xml';
[~,machineInit] = system('hostname');
machine = lower(deblank(cell2mat(regexp(machineInit, '^[^\.]+', 'match'))));

nevFilename = fullfile(dataPath, nevFilebase);
[~, trainParams, ~, ~] = readExperiment(xmlFile,subject,machine);

% In this context, nevFileBase matches nevFilesForTrain 
decoderFileLocationAndName = calibrateRewardAxisDecoderForBci('', nevFilename, {nevFilename}, trainParams, subject, true);
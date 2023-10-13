function [decoderFileLocationAndName] = trainDecoderFromNev(nevFilebase, xmlFile, subject)
% function [M0, M1, M2, channelsKeep, A, Q, C, R, beta, K] = trainKalmanDecoderFromNev(nev, trainParams, params, codes, datBase, nevBase,...
%         netLabels, gamma, channelNumbersUse, Qvalue, numLatents, velocityToCalibrateWith, includeBaseForTrain)

exGlobals;
dataPath = 'E:\';

[~,machineInit] = system('hostname');
machine = lower(deblank(cell2mat(regexp(machineInit, '^[^\.]+', 'match'))));

nevFilename = fullfile(dataPath, nevFilebase);
[~, trainParams, ~, ~] = readExperiment(xmlFile,subject,machine);

% In this context, nevFileBase matches nevFilesForTrain 
%decoderFileLocationAndName = calibrateMultipleAxisDecoderForBci('', nevFilename, {nevFilename}, trainParams, subject, true);
% Should be changed to match the file used with the nevFileBase
nevFilesForTrain = {
    'E:\satchel\Sa230831_s734_kalmanVelocityBci_0002.nev',...
    'E:\satchel\Sa230831_s734_kalmanVelocityBci_0003.nev',...
    'E:\satchel\Sa230831_s734_kalmanVelocityBci_0004.nev',...
    'E:\satchel\Sa230831_s734_kalmanVelocityBci_0005.nev',...
    'E:\satchel\Sa230831_s734_kalmanVelocityBci_0006.nev'
};
disp('whoops')
decoderFileLocationAndName = calibrateKalmanVelocityAndIntuitiveAxisAnnulusDecoderForBci('', nevFilename, nevFilesForTrain, trainParams, subject, true);

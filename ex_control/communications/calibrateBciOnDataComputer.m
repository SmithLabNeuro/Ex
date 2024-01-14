function [decoderTrainedOut, didTrainDecoder] = calibrateBciOnDataComputer(e, bciSocket)

global params socketsDatComp codes outfilename allCodes;
persistent blockLast decoderTrained

didTrainDecoder = false;

% setup block info
lenAC = length(allCodes);
thisBlock = e(1).currentBlock;

if e(1).trialCounter == 1 && thisBlock == 1
    blockLast = thisBlock; % reset here
    if strcmp(e(1).bciCalibration_bciDecoderFile, 'trainNew')
        % only reset if you're training a new decoder
        decoderTrained = false(size(e(1).bciCalibration_trainAfterBlock));
    elseif lenAC < 2 % check if even a trial has run
        % (or actually reset if the old one hasn't been loaded yet
        decoderTrained = false(size(e(1).bciCalibration_trainAfterBlock));
    end
end
if strcmp(e(1).bciCalibration_bciDecoderFile, 'trainNew')

    % checks if training should happen after last block (assuming we're at a
    % block change)
    segmentParamTraining = blockLast == e(1).bciCalibration_trainAfterBlock;
    trainAfterLastBlock = any(segmentParamTraining);
    % checks if the training has already happened
    lastBlockDecoderIsTrained = trainAfterLastBlock && decoderTrained(segmentParamTraining);
    
    if blockLast ~= thisBlock
        % flag for a block change, which is where training happens
        blockChange = true;
        blockLast = thisBlock;
    else
        blockChange = false;
    end
    if blockChange && trainAfterLastBlock && ~lastBlockDecoderIsTrained
        % we're set to train cuz we have some data
        rcAck = sendMessageWaitAck(socketsDatComp, 'trainDecoder');
        %     relNevFilepath = receiveMessageSendAck(socketsDatComp);
        %     writeExperimentInfoToDatabase([], [], outfilename, 'neural_output_name', relNevFilepath)
        nevFilenames = readExperimentInfoFromDatabase(outfilename, 'neural_output_name');
        nevFilenamesSplit = strsplit(nevFilenames{1}, '\n');
        % the base need to be read in because it has all the parameters!
        nevBaseFilename = [nevFilenamesSplit{1}, '.nev'];
        nevFilenamesUse = nevFilenamesSplit(e(1).bciCalibration_trainNextBlockDecoderUsingSegments{segmentParamTraining});
        nevFilesForTraining = [strjoin(nevFilenamesUse, '.nev\n'), '.nev'];
        
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).decoderCalibrationFunctionName);
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).bciDecoderParamFile);
        
        % receive the BCI parameters to send on with the new NEV file (after
        % training has been completed)
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(paramStructChars, 'startSendingAsciiParameters')
            keyboard
        end
        bciParamStructAsEvalChar = [];
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        while ~strcmp(paramStructChars, 'endSendingAsciiParameters')
            bciParamStructAsEvalChar = [bciParamStructAsEvalChar, paramStructChars];
            paramStructChars = receiveMessageSendAck(socketsDatComp);
        end
        
        rcAck = sendMessageWaitAck(socketsDatComp, nevBaseFilename);
        rcAck = sendMessageWaitAck(socketsDatComp, nevFilesForTraining);
        
        maxTimeToTrain = 1000; % seconds
        try
            decoderFilenameForBciComputer = receiveMessageSendAck(socketsDatComp, maxTimeToTrain);
        catch err
            if strcmp(err.identifier, 'communication:waitForData:communicationTimeoutWithDataComputer')
                keyboard
                decoderFilenameForBciComputer = receiveMessageSendAck(socketsDatComp, maxTimeToTrain);
            else
                rethrow(err);
            end
        end
        
        recStarted = startNevRecording(socketsDatComp, e(1).xmlFile,outfilename);
        sendStruct(struct('bciDecodeParameterFile', decoderFilenameForBciComputer));
        bciParamStructEncodedForNev = double(bciParamStructAsEvalChar)+256;
        sendCode(bciParamStructEncodedForNev);
        notesUpdate = sprintf('decoder trained on segments/nev files [%s]: %s', num2str(e(1).bciCalibration_trainNextBlockDecoderUsingSegments{segmentParamTraining}), decoderFilenameForBciComputer);
        writeExperimentInfoToDatabase([], [], outfilename, 'extra_notes', decoderFilenameForBciComputer)
        
        if recStarted
            sendMessageWaitAck(bciSocket, 'decoderParameterFile');
            sendMessageWaitAck(bciSocket, decoderFilenameForBciComputer);
            decoderTrained(blockLast-1 == e(1).bciCalibration_trainAfterBlock) = true;
        else
            keyboard
        end
        didTrainDecoder = true;
%         sendCode(result);
    end
elseif strcmp(e(1).bciCalibration_bciDecoderFile, 'useLastTrained')
    if ~all(decoderTrained)
        % this assumes no more training is needed
        decoderTrained = true(size(e(1).bciCalibration_trainAfterBlock));
        decoderFilenameForBciComputer = '';
        lastExperimentId = outfilename;
        while isempty(decoderFilenameForBciComputer)
            currentExperimentId = readExperimentInfoFromDatabase(lastExperimentId, 'rowid');
            lastExperimentId = currentExperimentId{1}-1;
            previousExperimentDecoders = readExperimentInfoFromDatabase(lastExperimentId, 'extra_notes');
            previousExperimentDecodersByCell = strsplit(previousExperimentDecoders{1}, '\n');
            decoderFilenameForBciComputer = previousExperimentDecodersByCell{end};
        end
        sendStruct(struct('bciDecodeParameterFile', decoderFilenameForBciComputer));
        
        rcAck = sendMessageWaitAck(socketsDatComp, 'sendDecoderParameters');
        % receive the BCI parameters to send on with the new NEV file (after
        % training has been completed)
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).bciDecoderParamFile);
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(paramStructChars, 'startSendingAsciiParameters')
            keyboard
        end
        bciParamStructAsEvalChar = [];
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        while ~strcmp(paramStructChars, 'endSendingAsciiParameters')
            bciParamStructAsEvalChar = [bciParamStructAsEvalChar, paramStructChars];
            paramStructChars = receiveMessageSendAck(socketsDatComp);
        end
        
        % making sure things started recording
        recStatus = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(recStatus, 'recording')
            keyboard
            %             recStarted = true;
        else
            recStarted = true;
        end
        
        % send the parameters to the NEV
        bciParamStructEncodedForNev = double(bciParamStructAsEvalChar)+256;
        sendCode(bciParamStructEncodedForNev);
        writeExperimentInfoToDatabase([], [], outfilename, 'extra_notes', decoderFilenameForBciComputer)
        
        % send the decoder parameters to the BCI computer!
        fprintf('sending parameter file %s\n', decoderFilenameForBciComputer)
        if recStarted
            sendMessageWaitAck(bciSocket, 'decoderParameterFile');
            sendMessageWaitAck(bciSocket, decoderFilenameForBciComputer);
        else
            keyboard
        end
        
        % basically restart this trial
        didTrainDecoder = true;
%         sendCode(result);
%         return
    end
elseif strcmp(e(1).bciCalibration_bciDecoderFile, 'trainNewOnPreviousNevFiles')
    if ~all(decoderTrained)
        % this assumes no more data needed before training starts
        
        % we're set to train cuz we have some data
        rcAck = sendMessageWaitAck(socketsDatComp, 'trainDecoder');
        
        lastExperimentId = outfilename;
        nevFilenamesSplit = {''};
        while isempty(nevFilenamesSplit{1})
            currentExperimentId = readExperimentInfoFromDatabase(lastExperimentId, 'rowid');
            
            lastExperimentId = currentExperimentId{1}-1;
            previousExperimentNevFilesTaskAndDecoders = readExperimentInfoFromDatabase(lastExperimentId, 'neural_output_name', 'task', 'extra_notes');
            if ~strcmp(previousExperimentNevFilesTaskAndDecoders{2}, e(1).bciCalibration_bciTaskToTrainOn)
                continue
            end
            nevFilenamesSplit = strsplit(previousExperimentNevFilesTaskAndDecoders{1}, '\n');
            decoderNamesCell = strsplit(previousExperimentNevFilesTaskAndDecoders{3}, '\n');
            lastDecoderRelPath = decoderNamesCell{end};
        end
        
        % the base need to be read in because it has all the parameters!
        nevBaseFilename = [nevFilenamesSplit{1}, '.nev'];
        % train on everything for now
        nevFilenamesUse = nevFilenamesSplit;
        nevFilesForTraining = [strjoin(nevFilenamesUse, '.nev\n'), '.nev'];
        
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).decoderCalibrationFunctionName);
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).bciDecoderParamFile);
        
        % receive the BCI parameters to send on with the new NEV file (after
        % training has been completed)
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(paramStructChars, 'startSendingAsciiParameters')
            keyboard
        end
        bciParamStructAsEvalChar = [];
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        while ~strcmp(paramStructChars, 'endSendingAsciiParameters')
            bciParamStructAsEvalChar = [bciParamStructAsEvalChar, paramStructChars];
            paramStructChars = receiveMessageSendAck(socketsDatComp);
        end
        
        rcAck = sendMessageWaitAck(socketsDatComp, nevBaseFilename);
        rcAck = sendMessageWaitAck(socketsDatComp, nevFilesForTraining);
        
        
        maxTimeToTrain = 1000; % seconds
        
        try
            nextInfo = receiveMessageSendAck(socketsDatComp, maxTimeToTrain);
        catch err
            if strcmp(err.identifier, 'communication:waitForData:communicationTimeoutWithDataComputer')
                keyboard
                nextInfo = receiveMessageSendAck(socketsDatComp, maxTimeToTrain);
            else
                rethrow(err);
            end
        end
        
        % kinda hacky way to see if the calibration file is asking
        % something of us
        if strcmp(nextInfo, 'kalmanDecoderForPlane')
            sendMessageWaitAck(socketsDatComp, lastDecoderRelPath)
            decoderFilenameForBciComputer = receiveMessageSendAck(socketsDatComp, maxTimeToTrain);
        else
            decoderFilenameForBciComputer = nextInfo;
        end
        
        recStarted = startNevRecording(socketsDatComp, e(1).xmlFile,outfilename);
        sendStruct(struct('bciDecodeParameterFile', decoderFilenameForBciComputer));
        bciParamStructEncodedForNev = double(bciParamStructAsEvalChar)+256;
        sendCode(bciParamStructEncodedForNev);
%         notesUpdate = sprintf('decoder trained on segments/nev files [%s]: %s', num2str(e(1).bciCalibration_trainNextBlockDecoderUsingSegments{segmentParamTraining}), decoderFilenameForBciComputer);
        writeExperimentInfoToDatabase([], [], outfilename, 'extra_notes', decoderFilenameForBciComputer)
        
        if recStarted
            sendMessageWaitAck(bciSocket, 'decoderParameterFile');
            sendMessageWaitAck(bciSocket, decoderFilenameForBciComputer);
            %             end
            % prevent more training from happening
            decoderTrained = true(size(e(1).bciCalibration_trainAfterBlock));
        else
            keyboard
        end
        didTrainDecoder = true;
%         sendCode(result);
    end
elseif ischar(e(1).bciCalibration_bciDecoderFile)
    if ~contains(e(1).bciCalibration_bciDecoderFile, params.SubjectID)
        warning('Decoder "%s" might not have been trained on subject "%s"', e(1).bciCalibration_bciDecoderFile, params.SubjectID);
    end
    % if we already initialized stuff, then decoderTrained gets set to all
    % true, so we can skip this step and keep going
    if ~all(decoderTrained)
        % this assumes no more training is needed
        decoderTrained = true(size(e(1).bciCalibration_trainAfterBlock));
        decoderFilenameForBciComputer = e(1).bciCalibration_bciDecoderFile;
        
        sendStruct(struct('bciDecodeParameterFile', decoderFilenameForBciComputer));
        
        
        rcAck = sendMessageWaitAck(socketsDatComp, 'sendDecoderParameters');
        % receive the BCI parameters to send on with the new NEV file (after
        % training has been completed)
        rcAck = sendMessageWaitAck(socketsDatComp, e(1).bciDecoderParamFile);
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(paramStructChars, 'startSendingAsciiParameters')
            keyboard
        end
        bciParamStructAsEvalChar = [];
        paramStructChars = receiveMessageSendAck(socketsDatComp);
        while ~strcmp(paramStructChars, 'endSendingAsciiParameters')
            bciParamStructAsEvalChar = [bciParamStructAsEvalChar, paramStructChars];
            paramStructChars = receiveMessageSendAck(socketsDatComp);
        end
        
        % making sure things started recording
        recStatus = receiveMessageSendAck(socketsDatComp);
        if ~strcmp(recStatus, 'recording')
            keyboard
            %             recStarted = true;
        else
            recStarted = true;
        end
        
        % send the parameters to the NEV
        bciParamStructEncodedForNev = double(bciParamStructAsEvalChar)+256;
        sendCode(bciParamStructEncodedForNev);
        writeExperimentInfoToDatabase([], [], outfilename, 'extra_notes', decoderFilenameForBciComputer)
        
        % send the decoder parameters to the BCI computer!
        fprintf('sending parameter file %s\n', decoderFilenameForBciComputer)
        if recStarted
            sendMessageWaitAck(bciSocket, 'decoderParameterFile');
            sendMessageWaitAck(bciSocket, decoderFilenameForBciComputer);
        else
            keyboard
        end
        
        % basically restart this trial
        didTrainDecoder = true;
%         sendCode(result);
%         return
    end
else
    availableOptionsFor_bciCalibration_bciDecoderFile = sprintf(['bciCalibration_bciDecoderFile parameter can be:\n'...
        '- ''trainNew'' to go through the calibration process using NEVs currently acquired,\n'...
        '- the pathway (as a string) to a decoder as saved in the database to use that decoder,\n'...
        '- ''useLastTrained'' to grab the last trained decoder from the database (though use this with caution),\n'...
        '- ''trainNewOnPreviousNevFiles'' to go through the calibration process but using previously acquired NEVs']);
    error('''%s'' not implemented yet as an option for bciCalibration_bciDecoderFile parameter. Here''s what''s available:\n\n%s', e(1).bciCalibration_bciDecoderFile, availableOptionsFor_bciCalibration_bciDecoderFile)
end

decoderTrainedOut = decoderTrained;
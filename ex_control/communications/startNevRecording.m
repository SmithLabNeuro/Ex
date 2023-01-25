function [recStarted, recFileName] = startNevRecording(socketsDatComp, xmlFile, outfilename)
% does the communication dance with the data computer to record the
% experiment to an NEV
global params notes sqlDb
timeout = 10;
recStarted = false;

rc = sendMessageWaitAck(socketsDatComp, 'record', timeout);
if isempty(rc)
    return
end
rc = sendMessageWaitAck(socketsDatComp, params.SubjectID, timeout);
if isempty(rc)
    return
end

sessionNum = params.sessionNumber;
if ~isnan(sessionNum)
    rc = sendMessageWaitAck(socketsDatComp, num2str(sessionNum), timeout);
    if isempty(rc)
        return
    end
else
    return
end

[~,xmlBase,~] = fileparts(xmlFile);
rc = sendMessageWaitAck(socketsDatComp, xmlBase, timeout);
if isempty(rc)
    return
end

neuralOutName = receiveMessageSendAck(socketsDatComp, timeout);

try
    recordingCheck = receiveMessageSendAck(socketsDatComp, timeout);
    if ~strcmp(recordingCheck, 'recording')
        keyboard
    end
catch err
    if strcmp(err.identifier, 'communication:waitForData:communicationTimeoutWithDataComputer')
        keyboard
    else
        rethrow(err);
    end
end

if ~isempty(sqlDb)
    writeExperimentInfoToDatabase([], [], outfilename, 'neural_output_name', neuralOutName)

    if strcmp(notes(end), newline)
        notes = sprintf('%s%s\n', notes, neuralOutName);
    else
        notes = sprintf('%s\n%s\n', notes, neuralOutName);
    end
    sqlDb.exec(sprintf('UPDATE experiment_session SET notes = "%s" WHERE session_number = %d AND animal = "%s"', notes, params.sessionNumber, params.SubjectID));
end

% if it got here we assume the recording started
recStarted = true;
function [modelParams, updatedReturn, taskParamReceived] = processBciControlMessage(controlCompSocket, ctrlMsg, modelParams)

global params

updatedReturn = []; 
taskParamReceived = false;
if ~isempty(ctrlMsg)
    % check if there's a new decoder or if there's an adjusted
    % param! (Important for calibration steps!)
    if strcmp(ctrlMsg, 'decoderParameterFile')
        decoderParameterLocation = params.bciDecoderBasePathBciComputer;
        decoderParameterFileRelativePath = receiveMessageSendAck(controlCompSocket);
        decoderParameterFileRelativePath(decoderParameterFileRelativePath=='\') = '/';
        decoderParameterFileFullPath = fullfile(decoderParameterLocation, decoderParameterFileRelativePath);
        modelParams = load(decoderParameterFileFullPath);
        fprintf('loaded new parameters from %s\n', decoderParameterFileRelativePath)
    elseif strcmp(ctrlMsg, 'requestParameters')
        parameterName = receiveMessageSendAck(controlCompSocket);
        parameterValue = modelParams.(parameterName);
        sendMessageWaitAck(controlCompSocket, typecast(parameterValue,'uint8'));
    elseif strcmp(ctrlMsg, 'taskParameters')
        % Assumes task parameters will be sent as char arrays that will
        % then be converted into uint8 and then doubles.
        receivedParameters = receiveMessageSendAck(controlCompSocket);
        updatedReturn = double(uint8(receivedParameters));
        sprintf('got the parameters %s', sprintf('%i ', updatedReturn))
        sendMessageWaitAck(controlCompSocket, uint8(receivedParameters));
        taskParamReceived = true;
    elseif ~ischar(ctrlMsg)
        updatedReturn = typecast(uint8(ctrlMsg), 'double')';
        
    else
        updatedReturn = uint8(ctrlMsg);
    end
end
function [modelParams, updatedReturn] = processBciControlMessage(controlCompSocket, ctrlMsg, modelParams)

global params

updatedReturn = [];
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
    else
        %updatedReturn = typecast(uint8(ctrlMsg), 'double')';
        %updatedReturn = char(ctrlMsg);
        %updatedReturn = ctrlMsg;
        updatedReturn = typecast(uint8(ctrlMsg), 'double')';
    end
end
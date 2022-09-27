function [bciParamFile, subjectId] = receiveBciDecoderInformation(controlCompSocket)
% this function is responsible for receiving information about the BCI
% parameters

bciParamFile = receiveMessageSendAck(controlCompSocket);
subjectId = receiveMessageSendAck(controlCompSocket);

end


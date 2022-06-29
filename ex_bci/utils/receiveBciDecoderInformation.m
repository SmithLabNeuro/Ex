function [bciParamFile, subjectId] = receiveBciDecoderInformation(controlCompSocket)
% this function is responsible for receiving information about the BCI
% parameters
% bciReadyMsg = checkForControlMsgs(controlCompSocket);
% if strcmp(bciReadyMsg, 'readyBciParamFile')
sockets.receiver = controlCompSocket;
sockets.sender = controlCompSocket;
    bciParamFile = receiveMessageSendAck(sockets);
    subjectId = receiveMessageSendAck(sockets);
% else
%     bciParamFile = [];
%     subjectId = [];
% end

end


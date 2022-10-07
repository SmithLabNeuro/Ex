function controlCompSocket = initializeControlComputerCommunication()

%% function initializes communication with control computer

global params
exGlobals;
matlabUDP2('all_close');
controlCompSocketSendReceive = matlabUDP2('open', params.bci2controlIP, params.control2bciIP, params.control2bciSocket);
controlCompSocket.sender = controlCompSocketSendReceive;
controlCompSocket.receiver = controlCompSocketSendReceive;

controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
while ~strcmp(controlMsg, 'ready')
    pause(0.1)
    controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
end
matlabUDP2('send', controlCompSocket.sender, 'prepared');
disp('prepared')
% flush the buffer in case multiple readies were sent
while ~isempty(controlMsg)
    pause(0.1);
    controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
end
% acknowledge flush
matlabUDP2('send',controlCompSocket.sender, 'flushed');
disp('flushed and set to go')
end
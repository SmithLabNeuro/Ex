function controlCompSocket = initializeControlComputerCommunication()

%% function initializes communication with control computer

global params
exGlobals;
matlabUDP2('all_close');
controlCompSocket = matlabUDP2('open', params.bci2controlIP, params.control2bciIP, params.control2bciSocket);
controlMsg = matlabUDP2('receive',controlCompSocket);
while ~strcmp(controlMsg, 'ready')
    pause(0.1)
    controlMsg = matlabUDP2('receive',controlCompSocket);
end
matlabUDP2('send', controlCompSocket, 'prepared');
disp('prepared')
% flush the buffer in case multiple readies were sent
while ~isempty(controlMsg)
    pause(0.1);
    controlMsg = matlabUDP2('receive',controlCompSocket);
end
% acknowledge flush
matlabUDP2('send',controlCompSocket, 'flushed');
disp('flushed and set to go')
end
function [controlCompSocket, controlMsg] = checkIfControlComputerCommunicationNeedsReinitialization(controlCompSocket)
% checks if 'ready' is received, which would indicate a
% reinitialization, but pops out the controlMsg if it's not
% 'ready', so that downstream functions can use it. This is
% basically a spruced up 'receive' call to matlabUDP2.
controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
if strcmp(controlMsg, 'ready')
    controlCompSocket = initializeControlComputerCommunication();
elseif ~isempty(controlMsg)
    % this function sends an ack for any other message
    matlabUDP2('send', controlCompSocket.sender, 'ack');
end
end
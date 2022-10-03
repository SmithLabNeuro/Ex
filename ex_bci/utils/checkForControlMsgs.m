function controlMsg = checkForControlMsgs(controlCompSocket)
    % this is basically a 'receive' call to matlabUDP2, but with a nicer
    % function name
    controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
    if ~isempty(controlMsg)
        matlabUDP2('send', controlCompSocket.sender, 'ack');
    end
end
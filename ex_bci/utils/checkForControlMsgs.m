function controlMsg = checkForControlMsgs(controlCompSocket)
    % this is basically a 'receive' call to matlabUDP2, but with a nicer
    % function name
    ind = 0;
    controlMsg = [];
    while matlabUDP2('check', controlCompSocket.receiver)
        ind = ind+1;
        controlMsg = matlabUDP2('receive',controlCompSocket.receiver);
        if ~isempty(controlMsg)
            matlabUDP2('send', controlCompSocket.sender, 'ack');
        end
    end
    if ind>1
        disp(ind);
    end
end
function rcvd = sendMessageWaitAck(sockets, msgToSend, timeout)
% function sends a message over UDP to a socket and awaits an 'ack'
% response
if nargin<3
    timeout = 10; % seconds
end

if isunix
    matlabUDP2('send', sockets.sender, msgToSend)
    rcvd = waitForData(sockets.receiver, 'ack', timeout);
elseif ispc
    if ischar(msgToSend)
        msgToSend = uint8(msgToSend);
    end

    sockets.sender(msgToSend);
    msg = char(sockets.receiver());
    while isempty(msg)
        pause(0.1)
        msg = char(sockets.receiver())';
    end
    
    if ~strcmp(msg, 'ack')
        keyboard
    end
end


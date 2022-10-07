function rcvd = receiveMessageSendAck(sockets, timeout)
% function receives a message over UDP on a socket and sends an 'ack' in
% response
if nargin<2
    timeout = 10; % seconds
end

if isunix
    rcvd = waitForData(sockets.receiver, '', timeout);
    matlabUDP2('send', sockets.sender, 'ack');
elseif ispc
    rcvd = char(sockets.receiver())';
    while isempty(rcvd)
        pause(0.1)
        rcvd = char(sockets.receiver())';
        if length(rcvd) == sockets.receiver.MaximumMessageLength
            warning('a message longer than the maximum buffer length may have been received; if so, this can be changed for sockets.receiver with the MaximumMessageLength parameter so there''s hope');
        end
    end
    %     disp(msg)
    sockets.sender(uint8('ack'));
end

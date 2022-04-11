function rcvd = sendMessageWaitAck(sockets, msgToSend, timeout)
% function sends a message over UDP to a socket and awaits an 'ack'
% response
if nargin<3
    timeout = 10; % seconds
end

matlabUDP2('send', sockets.sender, msgToSend)
rcvd = waitForData(sockets.receiver, 'ack', timeout);


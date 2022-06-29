function rcvd = receiveMessageSendAck(sockets, timeout)
% function receives a message over UDP on a socket and sends an 'ack' in
% response
if nargin<2
    timeout = 10; % seconds
end

rcvd = waitForData(sockets.receiver, '', timeout);
matlabUDP2('send', sockets.sender, 'ack');

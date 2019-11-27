function msg(m,args)
% function msg(msg,args)
% ex helper function, sends a message to the display and then moves on
%
% msg: a string message
% args: if present, the msg variable can use % replacement like fprintf and
%   these are the arguments
%
% examples:
% > msg('obj_on 3');
% > msg('obj_on %i',3);

global debug sockets;

s = '';
if nargin > 1
    s = sprintf(m, args);
else
    s = sprintf(m);
end

if debug
    fprintf('Sent: %s\n', s);
end

matlabUDP2('send',sockets(1), s);


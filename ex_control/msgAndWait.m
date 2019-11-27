function msgAndWait(m,args,timeout)
% function msg(msg,args)
% ex helper function, sends a message to the display and then waits for a
% response from the display
%
% msg: a string message
% args: if present, the msg variable can use % replacement like fprintf and
%   these are the arguments
%
% examples:
% > msgAndWait('obj_on 3');
% > msgAndWait('obj_on %i',3);

if nargin > 1
    m = sprintf(m, args);
end

if nargin<3, timeout = 1; end; %changed default from 10 to 1 -acs14mar2016

msg(m);

first_arg = strtok(m);

start = tic;
rcvd = '';
while ~strcmpi(rcvd,first_arg)
    if toc(start)>timeout, %very unlikely error - would require the message buffer getting filled as fast as you empty it... with undesired messages...
        error('msgAndWait:aborted','Timeout waiting for specific return message');
    end;
    rcvd = waitFor(timeout); %changed from a hardcoded 1 to 'timeout' -acs14mar2016
end;

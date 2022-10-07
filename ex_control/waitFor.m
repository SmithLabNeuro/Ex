    function s = waitFor(timeout)
% function waitFor(timeout)
%
% just waits for any input from the display, and then returns. 
%
% Used by msgAndWait to detect return messages. -acs15mar2016

global debug params sockets;

if nargin<1, timeout = 10; end;

thisStart = tic;
warnOn = false;
while true
    loopTop = GetSecs;
    
    if matlabUDP2('check',sockets(1))
        s = matlabUDP2('receive',sockets(1));
        if debug
            fprintf('Rcvd: %s\n', s);
        end;
        if strcmp(s,'abort'), %added to potentially help with hangs. -ACS 03Sep2013
            error('waitFor:aborted','Abort signal received');
        end;
        break; %message received
    elseif toc(thisStart)>timeout
        error('waitFor:aborted','Timeout waiting for return message');
    else
        s = '';
    end;
    chk = GetSecs-loopTop;
    if (chk)>params.waitForTolerance, warning('waitFor:tooSlow','waitFor exceeded latency tolerance - %s - by %0.01f msecs',datestr(now), 1000*(chk - params.waitForTolerance)); warnOn = true;  end; %warn tolerance exceeded -acs22dec2012
end

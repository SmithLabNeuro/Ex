function rcvd = waitForData(socketRec, msgToReceive, timeout)
% function 
st = tic;
rcvd = '';

while toc(st) < timeout && ((~strcmp(rcvd, msgToReceive) && ~isempty(msgToReceive)) || (isempty(msgToReceive) && isempty(rcvd)))
    if matlabUDP2('check',socketRec)
        rcvd = matlabUDP2('receive',socketRec);
    end
    pause(0.001)
end

if ~isempty(msgToReceive) && ~strcmp(rcvd,msgToReceive)
    rcvd = '';
    errorSt = 'Communication with data computer failed--is it running recordex.m?';
    error('communication:waitForData:communicationFailWithDataComputer', errorSt)
elseif isempty(msgToReceive) && isempty(rcvd)
    % here if it reached a timeout
    rcvd = '';
    errorSt = 'expected to receive a response from the data computer that didn''t get received before the timeout...';
    error('communication:waitForData:communicationTimeoutWithDataComputer', errorSt)
end
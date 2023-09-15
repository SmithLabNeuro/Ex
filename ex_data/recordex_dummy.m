function recordex_dummy(connectionType)

global params codes
% do this before adding paths
exGlobals; % grab the global parameters

% addpath(genpath('C:\Users\smithlab\Ex\ex_control'));
% addpath(genpath('C:\Users\smithlab\spikesort\'))

if nargin<1
    connectionType = 'tcp';
end


% The first close (at least under UDP), seems to recognize a connection
% upon Matlab start (that existed from before Matlab started!) but not
% close it--this recognition makes the attempt at opening fail, because
% xippmex is in the wrong state to open new connections. As a result,
% recordex fails. In any case, a second 'close' seems to fix this...
% xippmex('close'); BUT! Things get more confusing, so see the comment in
% the try-catch below.
try
    % In a twist of confusion, when you assign an output to
    % xippmex('close'), xippmex actually errors if everything is fine
    % because apparently it doesn't assign an output when things are
    % correctly closed. If it *does* assign an output, apparently that's
    % an indication that things didn't close correctly, at which point we
    % need to close again. So, here, if an output is assigned it continues
    % to the second xippmex('close'), and if an output is not assigned but
    % we request it, it errors and as long as that output error was what
    % happened, recordex continues happily.
    statCheck = xippmex('close');
    pause(0.1); % pause needed I think to prevent crash
    xippmex('close');
catch error
    if ~contains(error.message, 'One or more output arguments not assigned during call to "xippmex".')
        rethrow(error)
    end
end

% go = waitForMessage;
status = false;
while ~status
     pause(2)
    disp('checking');
    try
         status = xippmex(connectionType);
    catch error
        if contains(error.message, sprintf('xippmex command string ''%s'' not recognized', connectionType))
            status = xippmex(connectionType);
            if ~status
                error('not sure what''s happening')
            end
        elseif ~contains(error.message, 'could not find NIP')
            rethrow(error)
        end
    end
end

disp('done checking');
% have to add an operator if this is a TCP connection
if strcmp(connectionType, 'tcp')
    operatorId = 129; % this it the last octet of the IP address of Ripple
    timeInLoop = 1;
    while true
        try
            xippmex('addoper', operatorId);
        catch err
            if isempty(strfind(err.message, 'operator not found'))
                rethrow(err)
            else
                % for some reason if it's running right as Trellis opens
                % this error appears, so we make sure to only show it if it
                % continually appears...
                if timeInLoop > 1 
                    fprintf('Trouble connecting to Trellis. Is it open?\n')
                end
                pause(0.2)
                timeInLoop = timeInLoop+1;
                continue % this keeps us in the loop
            end
        end
        % don't like this construction but it pops us out of the loop if
        % the connection succeeds
        fprintf('Connected to Trellis.\n')
        break 
    end
end

fprintf('xippmex running\n')

% open communication port to control computer; numbers are selected based
% on ports in exGlobals; at some point would be ideal to move these values
% into there
udpPortFromControlReceive = 4245;
udpPortFromControlSend = 4246;
controlIPAddress = '128.2.245.223';
maxNumBytes = 2000; % this value mimics matlabUDP2.h; should be enough
udpr = dsp.UDPReceiver('RemoteIPAddress', controlIPAddress, 'LocalIPPort', udpPortFromControlReceive, 'MaximumMessageLength', maxNumBytes);
udps = dsp.UDPSender('RemoteIPAddress', controlIPAddress, 'RemoteIPPort', udpPortFromControlSend);

socketsControlComm.receiver = udpr;
socketsControlComm.sender = udps;

defaultPauseTime = 1;
longPauseTime = 1;
pauseTime = defaultPauseTime;
dataPath = 'E:\';
bciDecodersBasePath = params.bciDecoderBasePathDataComputer;
bciDecodersParameterPath = fullfile(bciDecodersBasePath, params.bciDecoderXmlParamFolder);

disp('near the end');
%xippmex('close');
while true
    
     msg = receiveMessageSendAck(socketsControlComm);
    
    pause(1);
end
end

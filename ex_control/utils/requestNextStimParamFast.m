function [nextStimParam, overTime] = requestNextStimParamFast(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, bciSockets)

% Request and receive a uStim parameter using pre-uStim state from bci
% computer. This function is called during a brief computation period after
% pre-uStim period.
% During calibration trials, this function does not do anything.

if decoderTrained
    timePassed = 0;
    overTime = 0;
    maxWaitTime = 40; % ms

    % wait until getting a message sent from bci comp
    while ~matlabUDP2('check',bciSockets.receiver)
        pause(0.001) % 1ms
        timePassed = timePassed+1; % ms
        
        % if a message is not sent by bci comp within maxWaitTime, break
        if timePassed>maxWaitTime
            % nextStimParam = 1;
            nextStimParam = randsample(1:15504,1); % randomly choose one uStim pattern among all 5 elec uStim
            overTime = 1;
            disp('pre-uStim computation over time')
            break
        end
    end
    
    % update next uStim param with the received message if the message is
    % received within maxWaitTime
    if overTime ~= 1
        receivedMsg = '';
        % receive a message from a bci comp
        while strcmp(receivedMsg, 'ack')||isempty(receivedMsg)
            receivedMsg = matlabUDP2('receive',bciSockets.receiver);
        end
        receivedMsgCasted = typecast(uint8(receivedMsg),'double');
        
        % update  the uStim parameter value
        if bciAlgoType==4
            nextStimParam = receivedMsgCasted(1);
        elseif bciAlgoType==5
            nextStimParam = receivedMsgCasted(2);
        end
    end

    % flash the buffer in case multiple messages were sent
    flashBuffer(bciSockets);
end



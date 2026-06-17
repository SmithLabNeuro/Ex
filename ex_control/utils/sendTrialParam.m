function receivedMsgCasted = sendTrialParam(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, numValidPatterns, bciSockets)

if decoderTrained % i.e., closed-loop trials
    
    % wait until getting a message sent from bci comp
    timePassed = 0; % keep track of the time to take to get the response from bci comp
    while ~matlabUDP2('check',bciSockets.receiver)
        pause(0.01) % 10ms
        timePassed = timePassed+10; % ms
        if timePassed>100
            break
        end
    end

%     % receive the message
    if timePassed>100
        receivedMsgCasted = [-1 -1];
    else
        receivedMsg = '';
        while strcmp(receivedMsg, 'ack')||isempty(receivedMsg)
            receivedMsg = matlabUDP2('receive',bciSockets.receiver);
        end
        receivedMsgCasted = typecast(uint8(receivedMsg),'double');
    end
    % flash the buffer in case multiple messages were sent
    flashBuffer(bciSockets);
end



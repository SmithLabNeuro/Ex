function nextStimParam = requestNextStimParam(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, numValidPatterns, bciSockets)

% Request and receive a uStim parameter for the next trial from bci
% computer. This function is called at the end of a trial (which includes 
% both correct and incorrect trials as long as they are completed until uStim offset). 
% During calibration trials, this function does not do anything.

if decoderTrained % i.e., closed-loop trials
    
    % wait until getting a message sent from bci comp
    timePassed = 0; % keep track of the time to take to get the response from bci comp
    while ~matlabUDP2('check',bciSockets.receiver)
        pause(0.01) % 10ms
        timePassed = timePassed+10; % ms
    end

    % receive the message
    receivedMsg = '';
    while strcmp(receivedMsg, 'ack')||isempty(receivedMsg)
        receivedMsg = matlabUDP2('receive',bciSockets.receiver);
    end
    receivedMsgCasted = typecast(uint8(receivedMsg),'double');

    % check some conditions to decide how to process the message
    if length(receivedMsgCasted)<2 % not enough length of Msg to assign nextStimParam
        isGreedy = 2; % invalid trial
    else
        isGreedy = receivedMsgCasted(2);
    end
    if receivedMsgCasted(1)>numValidPatterns(bciAlgoType) % selected pattern is out of the valid pattern index range
        isGreedy = 2; % invalid trial
    end
    if length(receivedMsgCasted)<3 % not enough length of Msg to assign bciTrialCnt
        bciTrialCnt = 0; % invalid trial
        isGreedy = 2; % invalid trial
    else
        bciTrialCnt = receivedMsgCasted(3);
    end 
    
    if isGreedy == 2 % trial where some communication error happened. so don't update next uStim parameter value
        rms(1) = receivedMsgCasted(1);
    else
        if bciAlgoType==1
            nextStimParam(1) = receivedMsgCasted(1);
        elseif bciAlgoType==2
            nextStimParam(2) = receivedMsgCasted(1);
        elseif bciAlgoType==3
            nextStimParam(3) = receivedMsgCasted(1);
        elseif bciAlgoType==4
            nextStimParam(4) = receivedMsgCasted(1);
        elseif bciAlgoType==5
            nextStimParam(5) = receivedMsgCasted(1);
        end
    end
    
    % save the closed-loop update related variables
    bciParams = struct();
    bciParams.nextStimParam = nextStimParam;
    bciParams.bciCalled = true;
    bciParams.isGreedy = isGreedy;
    bciParams.bciTrialCnt = bciTrialCnt;
    sendStruct(bciParams);
    
    % flash the buffer in case multiple messages were sent
    flashBuffer(bciSockets);
end



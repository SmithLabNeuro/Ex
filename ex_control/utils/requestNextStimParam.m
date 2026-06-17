function nextStimParam = requestNextStimParam(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, numValidPatterns, bciSockets, targetStateIndex)

% Request and receive a uStim parameter for the next trial from bci
% computer. This function is called at the end of a trial (which includes
% both correct and incorrect trials as long as they are completed until uStim offset).
% During calibration trials, this function does not do anything.
%
% For algo 3 (MiSO), updates the specific target's stim pattern using targetStateIndex.
% For other algos (1,2,4,5,6), updates all targets with the same value.

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
        % For algo 3, 7, and 8, update only the specific target's entry
        % For other algos, update all targets with the same value
        if bciAlgoType == 3
            % Update only the specific target for algo 3
            nextStimParam{targetStateIndex, 3} = receivedMsgCasted(1);
        elseif bciAlgoType == 7
            % Update only the specific target for algo 7 (per-target OMiSO+)
            nextStimParam{targetStateIndex, 7} = receivedMsgCasted(1);
        elseif bciAlgoType == 8
            % Update only the specific target for algo 8 (Greedy from Initial, per-target)
            nextStimParam{targetStateIndex, 8} = receivedMsgCasted(1);
        else
            % For algos 1,2,4,5,6: update all targets with the same value
            numTargets = size(nextStimParam, 1);
            for targetIdx = 1:numTargets
                nextStimParam{targetIdx, bciAlgoType} = receivedMsgCasted(1);
            end
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



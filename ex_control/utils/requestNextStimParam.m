function nextStimParam = requestNextStimParam(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, numValidPatterns, bciSockets)

if decoderTrained
    timePassed = 0;
    while ~matlabUDP2('check',bciSockets.receiver)
        pause(0.01) % 10ms
        timePassed = timePassed+10; % ms
    end

    receivedMsg = '';
    while strcmp(receivedMsg, 'ack')||isempty(receivedMsg)
        % receivedMsg = receiveMessageSendAck(bciSockets);
        receivedMsg = matlabUDP2('receive',bciSockets.receiver);
    end
    % receivedMsgCasted = typecast(receivedMsg,'double');
    receivedMsgCasted = typecast(uint8(receivedMsg),'double');
    % targetLatentVal = receivedMsgCasted(1);

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
    else
        bciTrialCnt = receivedMsgCasted(3);
    end
    
    if isGreedy == 2 % trial where some communication error happened. so don't update epsilon greedy values
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
        end
    end
    % nextStimParam = receivedMsgCasted(1);
    % targetLatentVal = receivedMsgCasted(3);
    % zPostMean = receivedMsgCasted(4);
    % bciTrialCnt = receivedMsgCasted(5);
    %disp(nextStimParam)
    %disp(isGreedy)
    
    % save the bci related variables
    bciParams = struct();
    bciParams.nextStimParam = nextStimParam;
    bciParams.bciCalled = true;
    bciParams.isGreedy = isGreedy;
    %bciParams.targetLatentVal = targetLatentVal;
    %bciParams.zPostMean = zPostMean;
    bciParams.bciTrialCnt = bciTrialCnt;
    sendStruct(bciParams);
    
    % flash the buffer in case multiple messages were sent
    flashBuffer(bciSockets);
%     receivedMsg = matlabUDP2('receive',bciSockets.receiver);
%     while ~isempty(receivedMsg)
%         pause(0.01) % 10ms
%         timePassed = timePassed+10; % ms
%         receivedMsg = matlabUDP2('receive',bciSockets.receiver);
%     end
% else
    % nextStimParam = nextStimParamOrig;
%     if bciAlgoType==1
%         nextStimParam(1) = nextStimParamOrig;
%     elseif bciAlgoType==2
%         nextStimParam(2) = nextStimParamOrig;
%     elseif bciAlgoType==3
%         nextStimParam(3) = nextStimParamOrig;
%     elseif bciAlgoType==4
%         nextStimParam(4) = nextStimParamOrig;
%     end
end



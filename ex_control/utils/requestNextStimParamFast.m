function [nextStimParam, overTime] = requestNextStimParamFast(decoderTrained, nextStimParamOrig, nextStimParam, bciAlgoType, bciSockets)

if decoderTrained
    % tic
    timePassed = 0;
    overTime = 0;
    while ~matlabUDP2('check',bciSockets.receiver)
        pause(0.001) % 1ms
        timePassed = timePassed+1; % ms
        if timePassed>40
            nextStimParam = 1;
            overTime = 1;
            disp('over 40ms')
            break
        end
    end
    
    if overTime ~= 1
        receivedMsg = '';
        while strcmp(receivedMsg, 'ack')||isempty(receivedMsg)
            receivedMsg = matlabUDP2('receive',bciSockets.receiver);
        end
        % receivedMsgCasted = typecast(receivedMsg,'double');
        receivedMsgCasted = typecast(uint8(receivedMsg),'double');
        if bciAlgoType==4
            nextStimParam = receivedMsgCasted(1);
        elseif bciAlgoType==5
            nextStimParam = receivedMsgCasted(2);
        end
        % targetLatentVal = receivedMsgCasted(1);
    end

    % flash the buffer in case multiple messages were sent
    flashBuffer(bciSockets);
%     receivedMsg = matlabUDP2('receive',bciSockets.receiver);
%     while ~isempty(receivedMsg)
%         receivedMsg = matlabUDP2('receive',bciSockets.receiver);
%     end
    % toc
end



function sendCode(code)
% function sendCode(code)
%
% Send a code (one or more integers) over the digital port using unixSendByte
% (linux), assuming the sendingCodes flag is set.
    
    global params thisTrialCodes trialTic
    
    if params.sendingCodes
        % new comedi-based Linux digital output function
        unixSendByte(code);
    end
    
    thisTrialCodes(end+1,:) = [code toc(trialTic)];
end

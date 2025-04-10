function sendCode(code)
% function sendCode(code)
%
% Send a code (one or more integers) over the digital port using unixSendByte
% (linux), assuming the sendingCodes flag is set.
    
    global params thisTrialCodes trialTic

    % some checks on codes that are sent
    assert(isnumeric(code),'sendCode: code must be numeric');
    assert(all(code>=0 & code<=2^16),'sendCode: code must be zero to 2^16');
    
    if params.sendingCodes
        % new comedi-based Linux digital output function
        unixSendByte(code);
    end
    
    thisTrialCodes(end+1:end+length(code),:) = [code(:) ones(length(code),1).*toc(trialTic)];
end

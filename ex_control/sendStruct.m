function sendStruct(s)
%function sendStruct(s)
%
% Send a struct as ascii over the digital port using unixSendByte
% (linux). This function does not support windows
%

global params thisTrialCodes trialTic

fields = fieldnames(s);

for i = 1:length(fields)
    val = s.(fields{i});
    val = num2str(val(:)'); % needs to be a row for the line below
    m = double([fields{i} '=' val ';'])+256;
    
    if params.sendingCodes
        % new comedi-based Linux digital output function
        unixSendByte(m);
    end
    
    thisTrialCodes(end+1:end+length(m),:) = [m(:) ones(length(m),1).*toc(trialTic)];
end


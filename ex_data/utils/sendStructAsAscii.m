function sendStructAsAscii(trainParams, sockets)

%function sendStructAsAscii(s)
%
% Send a struct as ascii over udp
%

fields = fieldnames(trainParams);

for i = 1:length(fields)
    val = trainParams.(fields{i});
    if iscell(val) % works for single layer cells, not nested ones
        for clInd = 1:length(val)
            valC = num2str(val{clInd}(:)'); % needs to be a row for the line below
            m = [fields{i} '{' num2str(clInd) '}=' valC ';'];
            
            sendMessageWaitAck(sockets, m);
        end
    else
        val = num2str(val(:)'); % needs to be a row for the line below
        m = [fields{i} '=' val ';'];
        
        sendMessageWaitAck(sockets, m);
    end
    

end

end
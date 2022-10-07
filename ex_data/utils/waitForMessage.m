function msg = waitForMessage(udpr, udps)
    msg = char(udpr())';
    while isempty(msg)
        pause(0.1)
        msg = char(udpr())';
        if length(msg) == udpr.MaximumMessageLength
            warning('a message longer than the maximum buffer length may have been received; if so, this can be changed for udpr with the MaximumMessageLength parameter so there''s hope');
        end
    end
%     disp(msg)
    udps(uint8('ack'));
end

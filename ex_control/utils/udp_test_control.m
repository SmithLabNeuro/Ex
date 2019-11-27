% udp_test_control.m
%
% Run the display side first, then this side. This program repeatedly sends
% messages from the control to the display to test UDP communication.

ipAddress = '192.168.1.10';
remotePort = 8866;
localPort = 8844;

delete(instrfind);
u = udp(ipAddress,'RemotePort',remotePort,'LocalPort',localPort);
fopen(u);

msg = 1;

while 1 
    fprintf(u,num2str(msg));

    while 1
        if get(u,'BytesAvailable') > 0            
            s = fgetl(u);
            if strcmp(s,num2str(msg))
                break;
            end
        end
    end
    
    disp(['Sent: ',num2str(msg)]);
    msg = msg + 1;
end

    
% udp_test_display.m
%
% Run this side first, then run the control side. This program just keeps
% listening for messages from the control side to test the UDP communication

ipAddress = '192.168.1.11';
remotePort = 8844;
localPort = 8866;

delete(instrfind);
u = udp(ipAddress,'RemotePort',remotePort,'LocalPort',localPort);
fopen(u);

display('Waiting for control side to test first message ...');

while 1
    if get(u,'BytesAvailable') > 0
        [s1 s] = strtok(fgetl(u));
        fprintf(u,s1);

        disp(['Received: ',s1]);
    end
end
%% test code for matlabUDP2, which allows UDP communication between an arbitrary number of computers
% Instructions:
% 1. start unix_matlabUDP2_test.m on the control computer
% 2. start unix_matlabUDP2_test_Display.m on the display computer
% 3. start unix_matlabUDP2_test_BCI.m on the BCI computer
% 4. Sit back and enjoy
%matlabUDP2('close',0);
matlabUDP2('all_close')

clear all
numiters = 100;

% IPs
ipA = '192.168.2.11'; % ip for control computer to BCI
ipB = '192.168.2.10'; % ip for BCI computer
ipC = '192.168.1.11'; % ip for control computer to display
ipD = '192.168.1.10'; % ip for display computer to control

% Ports
portA = 4595;
portB = 4596;
portC = 4597;
portD = 4598;
try
% open sockets
sockid = matlabUDP2('open',ipD,ipC,portC);

%% wait for control computer
donecode = [];
while ~isequal('readyDisplay',donecode)
    if matlabUDP2('check',sockid)
       donecode =  matlabUDP2('receive',sockid);   
    end
end
disp('Ready Display')
%% Test 2 Reliability and timing of communication with multiple computers

readvalsDisplay = zeros(numiters,1);
totaltimes = zeros(numiters,1);

for n = 1:numiters
    tic;
    while ~matlabUDP2('check',sockid)
    end
    temp = str2double(matlabUDP2('receive',sockid))
    readvalsDisplay(n) =    temp;
    matlabUDP2('send',sockid,num2str(readvalsDisplay(n)));
    totaltimes(n) = toc;
end


disp(['Mean Total Time: ', num2str(mean(totaltimes))])
disp(['Max Total Time: ', num2str(max(totaltimes))])
disp(['Min Total Time: ', num2str(min(totaltimes))])
 matlabUDP2('close',sockid);
catch
    matlabUDP2('close',sockid);
end
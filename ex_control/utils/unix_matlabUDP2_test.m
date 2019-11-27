%% test code for matlabUDP2, which allows UDP communication between an arbitrary number of computers
% Instructions:
% 1. start unix_matlabUDP2_test.m on the control computer
% 2. start unix_matlabUDP2_test_Display.m on the display computer
% 3. start unix_matlabUDP2_test_BCI.m on the BCI computer
% 4. Sit back and enjoy

%matlabUDP2('close',1)
%matlabUDP2('close',0)
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
%try
% open sockets

    sockid = matlabUDP2('open',ipA,ipB,portA);
    sockid2 = matlabUDP2('open',ipC,ipD,portC);

    %% Test 1 Reliability of packet receipt from BCI
    readycode = [];
    while ~isequal(readycode,'ready')
        if matlabUDP2('check',sockid)
           readycode =  matlabUDP2('receive',sockid);   
        end
    end
    disp('begin test')
    breakcode =num2str(numiters);
    testcode = '';
    timecheck = zeros(length(numiters),1);
    readvals = cell(length(numiters),1);
    readtimes = zeros(length(numiters),1);
    totaltime = zeros(length(numiters),1);
    count = 1;
    while ~isequal(testcode,breakcode)
        if matlabUDP2('check',sockid)
             
            tic; 
            matlabUDP2('check',sockid); 
            timecheck(count) =  toc;
             tstart =tic;
            testcode = matlabUDP2('receive',sockid);
           tstart2 = tic;
            readvals(count) =  {testcode};
            readtimes(count)=toc(tstart);  
            totaltime(count) =  toc(tstart2);
            count = count+1;
        end
    end
    % send done signal 2 BCI

    disp(['Mean Check Time: ', num2str(mean(timecheck))])
    disp(['Max Check Time: ', num2str(max(timecheck))])
    disp(['Min Check Time: ', num2str(min(timecheck))])
    disp(['Mean Read Time: ', num2str(mean(readtimes))])
    disp(['Max Read Time: ', num2str(max(readtimes))])
    disp(['Min Read Time: ', num2str(min(readtimes))])


%% Test 2 Reliability and timing of communication with multiple computers

matlabUDP2('send',sockid,'readyBCI');
matlabUDP2('send',sockid2,'readyDisplay');

readvalsBCI = cell(numiters,1);
totaltimes = zeros(numiters,1);
readvalsDisplay = cell(numiters,1);
readtimesDisplay = [];
testcodeBCI = [];
testcodeDisplay = [];
for n = 1:numiters
    tic;
    matlabUDP2('send',sockid,num2str(n));
    matlabUDP2('send',sockid2,num2str(n));
    while ~matlabUDP2('check',sockid)
    end
    readvalsBCI =  {matlabUDP2('receive',sockid)};
    while ~matlabUDP2('check',sockid2)
    end
    readvalsDisplay(n) = {matlabUDP2('receive',sockid2)};
    totaltimes(n) = toc;
end

disp(['Mean Total Time: ', num2str(mean(totaltimes))])
disp(['Max Total Time: ', num2str(max(totaltimes))])
disp(['Min Total Time: ', num2str(min(totaltimes))])
matlabUDP2('close',sockid);
matlabUDP2('close',sockid2);
%catch
%    matlabUDP2('close',sockid);
%    matlabUDP2('close',sockid2);
%end
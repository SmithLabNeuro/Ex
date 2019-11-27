% diode_loop_control.m
%
%

clear global params thisTrialCodes trialTic
global params thisTrialCodes trialTic

matlabUDP('close'); 
matlabUDP('open','192.168.1.11','192.168.1.10',4243);

exGlobals;

pause(1);

params.sendingCodes = 1;
thisTrialCodes = [];
trialTic = tic;

for I=1:1000
    disp(I)

    sendCode(codes.START_TRIAL);
    
    pause(0.1);
    
    sendCode(codes.STIM_ON);
    msgAndWait('diodeon');
    sendCode(codes.STIM_ON);
    
    pause(1);
    
    sendCode(codes.STIM_OFF);
    msgAndWait('diodeoff');
    sendCode(codes.STIM_OFF);
    
    pause(0.1);
    
    sendCode(codes.END_TRIAL);

    pause(0.1);

end

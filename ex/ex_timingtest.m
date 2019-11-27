function result = ex_timingtest(e)

global params codes behav;

phase = 0;
e.number = 20;
for i = 1:e.number
    rx = randi(600)-300;
    ry = randi(400)-200;
    msgAndWait('set %i grating %i %f %f %f %f %i %i %i %f',...
    [i 100  e.orientation phase e.spatial e.temporal e.centerx+rx e.centery+ry e.radius e.contrast]);
end
%msgAndWait('set 1 oval 0 %i %i %i %i %i %i %f',[e.fixX e.fixY e.fixRad 255 255 0 128]); %constant central fixation (yellow)

%msgAndWait('set 1 radialcheck 0 %i %i %i %f',[500 4 24 0.6]);
% msgAndWait('set 2 squarecheck 0 %i %i %i %i %i 255 %i %i %i 200 10 0 70 70',[30 8 255 255 255 0 0 0]);

msg('diode 1');

msgAndWait('ack');
msgAndWait(['obj_on ',num2str(1:e.number)]);
msgAndWait('timing_begin');
msgAndWait('diode_timing');
%msgAndWait('obj_on 2');
sendCode(codes.STIM_ON);

if ~waitForDisplay(e.fixX,e.fixY,params.fixWinRad)
    % failed to keep fixation
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);
    waitForMS(1000);
    result = codes.BROKE_FIX;
    msg('timing_end');
    frameDrop = str2double(waitFor());
    if frameDrop>0
        sendCode(codes.SHOWEX_TIMINGERROR);
        behav.frameDrop = [behav.frameDrop;frameDrop];     
    end
    return;
end

sendCode(codes.CORRECT);
msgAndWait('all_off'); % added MAS 2016/10/12 to turn off fix spot before reward
sendCode(codes.REWARD);
result = 1;

msg('timing_end');
frameDrop = str2double(waitFor());
if frameDrop>0
    sendCode(codes.SHOWEX_TIMINGERROR);
    result = codes.SHOWEX_TIMINGERROR;
    behav.frameDrop = [behav.frameDrop;frameDrop];
    fprintf('%d frame drops in the %d th blocks\n',frameDrop,e.currentBlock)
end
end


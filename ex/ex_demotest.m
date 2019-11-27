function result = ex_demotest(e)
% ex file: ex_driftchoice

global params codes behav allCodes;

phase = 0;
msgAndWait('set 2 grating %i %f %f %f %f %i %i %i %f',...
    [e.stimLength  e.orientation phase e.spatial e.temporal e.centerx e.centery e.radius e.contrast]);
msgAndWait('set 3 rect %i %i %i %i %i %i %i %i %i %i %i',...
    [e.stimLength  -e.centerx -e.centery e.width/2 e.height/2 0 200 200]);
msgAndWait('set 4 squarecheck 0 %i %i %i %i %i %i %i %i %i %i %i 0 %i %i -1',[e.squarer e.squaren e.squareColor1 e.squareColor2 e.seed e.stimx e.stimy]);
msgAndWait('set 5 radialcheck 0 %i %i %i %f %i %i %i',[e.screenYpix e.radialr e.radialn 255 -e.stimx e.stimy 0]);

msgAndWait('set 6 squarecheck %i %i %i %i %i %i %i %i %i %i %i %i 0 %i %i -1',[e.stimLength e.squarer e.squaren e.squareColor1 e.squareColor2 e.seed 0 0]);
msgAndWait('set 7 radialcheck 0 %i %i %i %f %i %i %i',[e.screenYpix e.radialr e.radialn e.radialAlpha 0 0 1]);

msgAndWait('set 1 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.fixRad 255 255 0]); %constant central fixation (yellow)
msg('diode 1');

msgAndWait('ack');
msgAndWait('obj_on 1');
sendCode(codes.FIX_ON);

if ~waitForFixation(1000,e.fixX,e.fixY,params.fixWinRad)
    % failed to achieve fixation
    sendCode(codes.IGNORED);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.IGNORED;
    return;
end

if ~waitForMS(e.preStimFix,e.fixX,e.fixY,params.fixWinRad)
    % hold fixation before stimulus comes on
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(1000);
    result = codes.BROKE_FIX;
    return;
end

msgAndWait('timing_begin');
msgAndWait('diode_timing');
msgAndWait('obj_on 2');
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
    return;
end
%msgAndWait('obj_off 2');
sendCode(codes.STIM_OFF);

msgAndWait('obj_on 3');
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
    return;
end
%msgAndWait('obj_off 3');
sendCode(codes.STIM_OFF);

msgAndWait('obj_on 4');
sendCode(codes.STIM_ON);

if ~waitForMS(500,e.fixX,e.fixY,params.fixWinRad)
    % failed to keep fixation
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);
    waitForMS(1000);
    result = codes.BROKE_FIX;
    msg('timing_end');
    return;
end
msgAndWait('obj_off 4');
sendCode(codes.STIM_OFF);

msgAndWait('obj_on 5');
sendCode(codes.STIM_ON);

if ~waitForMS(500,e.fixX,e.fixY,params.fixWinRad)
    failed to keep fixation
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);
    waitForMS(1000);
    result = codes.BROKE_FIX;
    msg('timing_end');
    return;
end
msgAndWait('obj_off 5');
sendCode(codes.STIM_OFF);

msgAndWait('obj_on 6 7');
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
    return;
end
sendCode(codes.STIM_OFF);


sendCode(codes.CORRECT);
msgAndWait('all_off'); % added MAS 2016/10/12 to turn off fix spot before reward
sendCode(codes.REWARD);

result = 1;

msg('timing_end');

end


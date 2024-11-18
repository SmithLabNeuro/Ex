function result = ex_kalmanVelocityDemoBci(e)

global params codes bciSockets exPrint;

[decoderTrained, didTrainDecoder] = calibrateBciOnDataComputer(e, bciSockets);
if didTrainDecoder
    result = codes.BACKGROUND_PROCESS_TRIAL;
    sendCode(result);
    return
end
if ~any(decoderTrained)
    assistIndToUse = 1;
else
    assistIndToUse = find(decoderTrained, 1, 'last') + 1;
end

if all(decoderTrained)
    movementTime = e(1).movementBciTime;
    holdHeMsDuringMovement = e(1).holdHeMsDuringMovement;
else
    movementTime = e(1).calibrationDecoderMovementTime;
    holdHeMsDuringMovement =  e(1).calibrationDecoderHoldHeMsDuringMovement;
end

centerToTargetBciScale = e(1).bciCalibration_centerToTargetBciScaleBySegment(assistIndToUse);
orthogonalBciScale = e(1).bciCalibration_orthogonalBciScaleBySegment(assistIndToUse);
centerToOutVelForAssist = e(1).bciCalibration_centerToOutVelForAssistBySegment(assistIndToUse);

angleTolerance = .5; 

% obj 1 is fix spot, obj 2 is stimulus, diode attached to obj 2, cursor is
% obj 3
fixObjId = 1;
targetObjectId = 2;
cursorObjId = 3;
msg(['diode ' num2str(targetObjectId)]);
msgAndWait('ack');

% Atari joystick requirements for holding down
joystickDevInds = connectToJoysticks(e(1).joystickName); % connecting to Atari joystick
if isfield(e(1), 'joystickXYDown')
    joystickXYDown = e(1).joystickXYDown;
else
    joystickXYDown = [-32768 0]; % push down for success
end
checkAll = true;
joystickAtariAcqPosArgs = {joystickXYDown(1), joystickXYDown(2), joystickDevInds, checkAll};

% fixation window and location
fixWinVals = params.fixWinRad;
fixX = e(1).fixX;
fixY = e(1).fixY;
fixationAcqArgs = {fixX, fixY, fixWinVals};

% Hall Effect joystick requirements for hold/twist
joystickHEXYHold = [0 0];
joystickHEAngHoldAtLeast = e(1).heJoystickMinHoldAngle;
distanceTolerance = e(1).joystickHoldTolerance;
joystickHEAcqPosArgs = {joystickHEXYHold(1), joystickHEXYHold(2), joystickHEAngHoldAtLeast, distanceTolerance, angleTolerance, e(1).bciMaxPixelDist};

% turn on fixation spot
msg('set %d oval 0 %i %i %i %i %i %i',[fixObjId, e(1).fixX e(1).fixY e(1).fixRad e(1).fixColor(1) e(1).fixColor(2) e(1).fixColor(3)]);
msgAndWait('obj_on %d', fixObjId);
sendCode(codes.FIX_ON);

% now wait for the subject to get into position (Atari, HE, fixation)
[conditional, perFuncConditional] = waitForEvent(e(1).timeToStaticJoystickGrab, {@joystickAtariAcquirePosition, @joystickHallEffectGrabCheck, @fixationAcquire}, {joystickAtariAcqPosArgs, joystickHEAcqPosArgs, fixationAcqArgs});

if ~conditional
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.IGNORED;
    sendCode(result);
    waitForMS(1000);
    return
end
% send code that he fixated
sendCode(codes.FIXATE);

% Monkey must keep holding down atari joystick and twist HE joystick for threeMechanicFixation milliseconds
joystickAtariHoldArgs = joystickAtariAcqPosArgs(1:3);
joystickHEHoldMsArgs = [joystickHEAcqPosArgs, cursorObjId, e(1).cursorRad, e(1).cursorColor, false, e(1).fixationJoystickHold];

fixationHoldArgs = fixationAcqArgs;

% now wait for subject to keep hold for threeMechanicFixation time
[overallSucc, perFuncConditional] = waitForEvent(e(1).fixationJoystickHold, {@joystickAtariHold, @joystickHallEffectHoldForMs, @fixationHold}, {joystickAtariHoldArgs, joystickHEHoldMsArgs, fixationHoldArgs});
if ~overallSucc
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.BROKE_FIX;
    sendCode(result);
    waitForMS(1000); % this is a timeout if he fails the holds
    return
end

% prep the cursor
cursorColorDisp = e(1).cursorColor;
startCursX = 0;
startCursY = 0;
cursorR = e(1).cursorRad;

% prep cursor for showing (NOTE: it's NOT on yet)
msg('set %d oval 0 %i %i %i %i %i %i', [cursorObjId startCursX startCursY cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

% choose a target location randomly around a circle
theta = deg2rad(e(1).targetAngle);
targetX = round(e(1).targetDistance * cos(theta)) + fixX;
targetY = round(e(1).targetDistance * sin(theta)) + fixY;
targWinCursRad = e(1).targWinCursRad;

% prep target for showing (NOTE: it's NOT on yet)
msg('set %d oval 0 %i %i %i %i %i %i',[targetObjectId targetX targetY e(1).targRad e(1).targColor(1) e(1).targColor(2) e(1).targColor(3)]);

% turn on target
msgAndWait('obj_on %d', targetObjectId);
sendCode(codes.TARG_ON);

% delayPeriodMs sets how long the target is on and the mechanics must be held before the cursor turns on allowing movement
joystickHEHoldMsArgs{end} = e(1).delayPeriodMs;
if ~waitForEvent(e(1).delayPeriodMs, {@joystickAtariHold, @joystickHallEffectHoldForMs, @fixationHold}, {joystickAtariHoldArgs, joystickHEHoldMsArgs, fixationHoldArgs})
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.FALSE_START;
    sendCode(result);
    waitForMS(1000); % this is a timeout if he fails the holds
    return
end

% turn off fixation point
msgAndWait('obj_off %d', fixObjId);
sendCode(codes.FIX_OFF);

% turn on cursor
msgAndWait('obj_on %d', cursorObjId);
sendCode(codes.CURSOR_ON)

% prep inputs for grabbing BCI-based-velocity
binSizeMs = 50;
pixelBoxLimit = e(1).bciMaxPixelDist;

% institute a freeze period
freezePeriod = e(1).freezePeriod;
joystickHEHoldMsArgs{end} = freezePeriod;
[heldForReaction, perFuncConditional, extraFuncOutputs] = waitForEvent(freezePeriod, {@joystickAtariHold, @joystickHallEffectHoldForMs}, {joystickAtariHoldArgs, joystickHEHoldMsArgs});

if ~heldForReaction
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.CURSOR_OFF);
    result = codes.FALSE_START;
    sendCode(result);
    waitForMS(1000); % this is a timeout if he reacts too fast
    return
end

% setup how subject can move the cursor to the target
initBciCursorPos = [0; 0];
initBciVelocity = [0; 0];
bciCursorMovementArgs = {bciSockets, binSizeMs, initBciCursorPos, initBciVelocity, pixelBoxLimit, targetX,targetY, e(1).targRad, cursorObjId, cursorR, cursorColorDisp, targWinCursRad, centerToTargetBciScale, orthogonalBciScale, centerToOutVelForAssist};
joystickHEHoldMsArgs{end} = holdHeMsDuringMovement;
heJoystickMinHoldAngleDuringMovement = e(1).heJoystickMinHoldAngleDuringMovement; % remove the need to twist during movement, but he can't knock the joystick!
joystickHEHoldMsArgs{3} = heJoystickMinHoldAngleDuringMovement; % reset the need to twist during the move
successUnlessFailBeforeHoldTime = true; % always a success unless he explicitly fails; and then always a success after the hold time
joystickHEHoldMsArgs{end+1} = successUnlessFailBeforeHoldTime;

% flush bci socket
ind = 0;
while matlabUDP2('check', bciSockets.sender)
    receivedMsg = matlabUDP2('receive', bciSockets.sender);
    ind = ind+1;
end

% Identify x-coord and y-coord for target opposite of intended
tmpStruct.decoderTrained = all(decoderTrained); % whether to use calib or bci params
sendStruct(tmpStruct);

[reachTarget, funcSuccesses] = waitForEvent(movementTime, {@bciMovementToTargetOptionalAssist, @joystickAtariHold, @joystickHallEffectHoldForMs}, {bciCursorMovementArgs, joystickAtariHoldArgs, joystickHEHoldMsArgs});
if ~reachTarget    
    if ~all(decoderTrained)
        result = codes.LATE_CHOICE;
    else
        result = codes.WRONG_TARG;
    end
    
    if funcSuccesses(1)<0
        % moved cursor too slow
        msgAndWait('all_off');
        sendCode(codes.CURSOR_OFF);
        result = codes.WRONG_TARG;
        sendCode(result);
        waitForMS(movementTime + 1000);  % timeout period (disincentivize him from skiqxpping)
    elseif funcSuccesses(3) < 0
        msgAndWait('all_off');
        sendCode(codes.CURSOR_OFF);
        result = codes.BROKE_FIX;
        sendCode(result);
        waitForMS(500);  % timeout period
    else
        % moved cursor too slow
        msgAndWait('all_off');
        sendCode(codes.CURSOR_OFF);
        sendCode(result);
        giveJuice(params.juiceX,params.juiceInterval,params.juiceTTLDuration*0.5);
        waitForMS(500);  % timeout period
    end
    return
end
sendCode(codes.CURSOR_OFF); % end of BCI period
sendCode(codes.ACQUIRE_TARG);
% CURSOR_OFF code needed for BCI computer to know the BCI is over

joystickHEAngHoldAtMost = e(1).heJoystickMaxHoldAngleForRelease;
joystickHEAcqPosArgs = {joystickHEXYHold(1), joystickHEXYHold(2), joystickHEAngHoldAtMost, distanceTolerance, angleTolerance, e(1).bciMaxPixelDist};

[reactedTooFast, perFuncConditionalFastReaction, extraFuncOutputsFastReaction] = waitForEvent(e(1).rtCatchMsForHeJoystickUntwist, {@joystickAtariAcquirePosition, @joystickHallEffectReleaseCheck}, {joystickAtariAcqPosArgs, joystickHEAcqPosArgs});
if e(1).rtCatchMsForHeJoystickUntwist && reactedTooFast
    % untwisted joystick too fast or let go of Atari
    msgAndWait('all_off');
    sendCode(codes.CURSOR_OFF);
    result = codes.BROKE_TARG;
    sendCode(result)
    waitForMS(500); % normal fail inter-trial-interval
    return
end

releaseTime = e(1).timeToStaticJoystickRelease - e(1).rtCatchMsForHeJoystickUntwist;
[conditional, perFuncConditional, extraFuncOutputs] = waitForEvent(releaseTime, {@joystickAtariAcquirePosition, @joystickHallEffectReleaseCheck}, {joystickAtariAcqPosArgs, joystickHEAcqPosArgs});
if ~e(1).timeToStaticJoystickRelease || conditional    
    sendCode(codes.REWARD);
    giveJuice(params.juiceX,params.juiceInterval,params.juiceTTLDuration);
    result = codes.CORRECT;
    sendCode(result);
else
    result = codes.NO_CHOICE;
    sendCode(result);
end
msgAndWait('all_off'); % turn off cursor/target after reward so he can see why he succeeded (or failed)
waitForMS(500); % intertrial interval
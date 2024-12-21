function result = ex_kalmanVelocityDemoBci(e)

global params codes bciSockets exPrint;

% This function checks weather or not to train a new decoder before
% starting this trial, and trains the decoder if nesecceray. 
% didTrainDecoder - boolean, will be true if the decoder was trained on
% this sepcific trial.
% decoderTrained - an array the size of e.bciCalibration_trainAfterBlock that tracks whether or
% not a decoder was arleady trained for a given block. 
[decoderTrained, didTrainDecoder] = calibrateBciOnDataComputer(e, bciSockets);
tmpStruct.decoderTrained = all(decoderTrained); % save decoder information to param file
sendStruct(tmpStruct);

% Set up trial params:
if didTrainDecoder
    % send code to record that a decoder has been trained
    result = codes.BACKGROUND_PROCESS_TRIAL; 
    sendCode(result);
    return
end

% assistIndToUse - this index defines the block, i.e, how much control vs 
% help the monkey gets for directing the cursor. 
if ~any(decoderTrained)
    assistIndToUse = 1;
else
    assistIndToUse = find(decoderTrained, 1, 'last') + 1;
end

% determines whether any calibration trials are left, 
% sets movement parameters to calibration or BCI values
if all(decoderTrained)
    movementTime = e(1).bciMovementTime;
    holdHeMsDuringMovement = e(1).holdHeMsDuringBciMovement;
else
    movementTime = e(1).calibrationMovementTime;
    holdHeMsDuringMovement =  e(1).holdHeMsDuringCalibrationMovement;
end

% How much control the monkey has is the target direction:
centerToTargetBciScale = e(1).bciCalibration_centerToTargetBciScaleBySegment(assistIndToUse);
% How much control the monkey has is the orthogonal to the target
% direction:
orthogonalBciScale = e(1).bciCalibration_orthogonalBciScaleBySegment(assistIndToUse);
% automonkey velocity
centerToOutVelForAssist = e(1).bciCalibration_centerToOutVelForAssistBySegment(assistIndToUse);
% When the cursor is angleTolerance of the way to the target, if it is not
% within the pie slice of the target the trial will fail. 
angleTolerance = e(1).angleTolerance; 

% showex objects: lowest object ID is drawn on top!
fixObjId = 1;
targetObjectId = 2;
cursorObjId = 3;
msg(['diode ' num2str(targetObjectId)]); % attach diode to target
msgAndWait('ack');

% Atari joystick requirements for holding down
joystickDevInds = connectToJoysticks(e(1).joystickName); % connecting to Atari joystick
joystickXYDown = e(1).joystickXYDown; % NaN ignores the Atari joystick, numbers requires a press
checkAll = true; % whether to check all dimensions of the Atari press
joystickAtariAcqPosArgs = {joystickXYDown(1), joystickXYDown(2), joystickDevInds, checkAll};

% fixation window and location
fixWinVals = params.fixWinRad;
fixX = e(1).fixX;
fixY = e(1).fixY;
fixationAcqArgs = {fixX, fixY, fixWinVals};
msg('set %d oval 0 %i %i %i %i %i %i',[fixObjId, e(1).fixX e(1).fixY e(1).fixRad e(1).fixColor(1) e(1).fixColor(2) e(1).fixColor(3)]);

% Hall Effect joystick requirements for hold/twist
joystickHEXYHold = [0 0];
joystickHEAngHoldAtLeast = e(1).heJoystickMinHoldAngle;
distanceTolerance = e(1).joystickHoldTolerance;
joystickHEAcqPosArgs = {joystickHEXYHold(1), joystickHEXYHold(2), joystickHEAngHoldAtLeast, distanceTolerance, angleTolerance, e(1).bciMaxPixelDist};

% Joystick holding args
joystickAtariHoldArgs = joystickAtariAcqPosArgs(1:3);
joystickHEHoldMsArgs = [joystickHEAcqPosArgs, cursorObjId, e(1).cursorRad, e(1).cursorColor, false, e(1).fixationJoystickHold];
fixationHoldArgs = fixationAcqArgs;

% prepare cursor arguments
cursorColorDisp = e(1).cursorColor;
startCursX = 0;
startCursY = 0;
cursorR = e(1).cursorRad;
% prep cursor for showing (NOTE: it's NOT on yet)
msg('set %d oval 0 %i %i %i %i %i %i', [cursorObjId startCursX startCursY cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

% prepare the target
theta = deg2rad(e(1).targetAngle);
targetX = round(e(1).targetDistance * cos(theta)) + fixX; % x position of target
targetY = round(e(1).targetDistance * sin(theta)) + fixY; % y position of target
targWinCursRad = e(1).targWinCursRad; % invisible window for successful target acquisition
% prep target for showing (NOTE: it's NOT on yet)
msg('set %d oval 0 %i %i %i %i %i %i',[targetObjectId targetX targetY e(1).targRad e(1).targColor(1) e(1).targColor(2) e(1).targColor(3)]);

% prepare BCI arguments
binSizeMs = e(1).bciBinsize;
pixelBoxLimit = e(1).bciMaxPixelDist;
freezePeriod = e(1).freezePeriod;
initBciCursorPos = [0; 0];
initBciVelocity = [0; 0];
bciCursorMovementArgs = {bciSockets, binSizeMs, initBciCursorPos, initBciVelocity, pixelBoxLimit, targetX,targetY, e(1).targRad, cursorObjId, cursorR, cursorColorDisp, targWinCursRad, centerToTargetBciScale, orthogonalBciScale, centerToOutVelForAssist};
% setup how subject can move the cursor to the target
successUnlessFailBeforeHoldTime = true; % always a success unless he explicitly fails; and then always a success after the hold time
joystickHEHoldDuringBci = [joystickHEAcqPosArgs, cursorObjId, e(1).cursorRad, e(1).cursorColor, false, holdHeMsDuringMovement, successUnlessFailBeforeHoldTime];
joystickHEHoldDuringBci{3} = e(1).heJoystickMinHoldAngleDuringMovement; % reset the need to twist during the move for BCI

%% TRIAL START

% Fixation aquisition: 
msgAndWait('obj_on %d', fixObjId);
sendCode(codes.FIX_ON);

% now wait for the subject to get into position (Atari, HE, fixation)
[acquiredFixation, ~] = waitForEvent(e(1).timeToStaticJoystickGrab, {@joystickAtariAcquirePosition, @joystickHallEffectGrabCheck, @fixationAcquire}, {joystickAtariAcqPosArgs, joystickHEAcqPosArgs, fixationAcqArgs});
if ~acquiredFixation
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.IGNORED;
    sendCode(result);
    waitForMS(e(1).punishmentTimeOutMs);
    return
end
% subject successfully fixated
sendCode(codes.FIXATE);

% Maintain fixation: 
% Monkey must keep holding down atari joystick and twist HE joystick for fixationJoystickHold milliseconds
[maintainedFixation, ~] = waitForEvent(e(1).fixationJoystickHold, {@joystickAtariHold, @joystickHallEffectHoldForMs, @fixationHold}, {joystickAtariHoldArgs, joystickHEHoldMsArgs, fixationHoldArgs});
if ~maintainedFixation
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.BROKE_FIX;
    sendCode(result);
    waitForMS(e(1).punishmentTimeOutMs); % this is a timeout if he fails the holds
    return
end

% turn on target
msgAndWait('obj_on %d', targetObjectId);
sendCode(codes.TARG_ON);

% delayPeriodMs sets how long the target is on and the mechanics must be held before the cursor turns on allowing movement
joystickHEHoldMsArgs{end} = e(1).delayPeriodMs; % update holdMS argument to proper time limit
if ~waitForEvent(e(1).delayPeriodMs, {@joystickAtariHold, @joystickHallEffectHoldForMs, @fixationHold}, {joystickAtariHoldArgs, joystickHEHoldMsArgs, fixationHoldArgs})
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.FALSE_START;
    sendCode(result);
    waitForMS(e(1).punishmentTimeOutMs); % this is a timeout if he fails the holds
    return
end

% turn off fixation point
msgAndWait('obj_off %d', fixObjId);
sendCode(codes.FIX_OFF);

% turn on cursor
msgAndWait('obj_on %d', cursorObjId);
sendCode(codes.CURSOR_ON)

% institute a freeze period - for proper BCI decoder initialization
joystickHEHoldMsArgs{end} = freezePeriod;
[heldForFreeze, ~, ~] = waitForEvent(freezePeriod, {@joystickAtariHold, @joystickHallEffectHoldForMs}, {joystickAtariHoldArgs, joystickHEHoldMsArgs});
if ~heldForFreeze
    % moved cursor too or broke fixation too early
    msgAndWait('all_off');
    sendCode(codes.CURSOR_OFF);
    result = codes.FALSE_START;
    sendCode(result);
    waitForMS(e(1).punishmentTimeOutMs); % this is a timeout if he reacts too fast
    return
end

%% BCI portion
% flush bci socket
ind = 0;
while matlabUDP2('check', bciSockets.sender)
    receivedMsg = matlabUDP2('receive', bciSockets.sender); % ???
    ind = ind+1;
end

% Now the cursor begins to move
[reachTarget, funcSuccesses] = waitForEvent(movementTime, {@bciMovementToTargetOptionalAssist, @joystickAtariHold, @joystickHallEffectHoldForMs}, {bciCursorMovementArgs, joystickAtariHoldArgs, joystickHEHoldDuringBci});
if ~reachTarget    
    if ~all(decoderTrained)
        result = codes.LATE_CHOICE; % calibration result
    else
        result = codes.WRONG_TARG; % BCI result 
    end
    
    if funcSuccesses(2)<0 ||  funcSuccesses(3)<0
        % aborted the task by releasing a joystick
        msgAndWait('all_off');
        sendCode(codes.CURSOR_OFF);
        result = codes.BROKE_TASK;
        sendCode(result);
        waitForMS(e(1).punishmentTimeOutMs);  % timeout period
    else
        % moved cursor too slow
        msgAndWait('all_off');
        sendCode(codes.CURSOR_OFF);
        sendCode(codes.REWARD);
        giveJuice(params.juiceX,params.juiceInterval,params.juiceTTLDuration*0.5); % partial reward for trying
        sendCode(result);
        waitForMS(e(1).interTrialInterval);
    end
    return
end
sendCode(codes.CURSOR_OFF); % end of BCI period, CURSOR_OFF code needed for BCI computer to know the BCI is over
sendCode(codes.ACQUIRE_TARG);

% if he made it to the target, he succeeds at the task!
sendCode(codes.REWARD);
giveJuice(params.juiceX,params.juiceInterval,params.juiceTTLDuration);
result = codes.CORRECT;
sendCode(result);
msgAndWait('all_off');
waitForMS(e(1).interTrialInterval); % intertrial interval
return
function success = joystickHallEffectShowCursor(loopStart, loopNow, cursorR, cursorColor, pixelDistForMaxJoystickPos, msHold)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes

pixBoxLimit = pixelDistForMaxJoystickPos;


% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

% compute the new cursor position (and don't let it get out of the bounding
% box)
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
cursorPos = [xVal*pixBoxLimit, yVal*pixBoxLimit]; 

posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);


% compute how close the cursor is to the target

if keyboardEvents()
    success = -1;
end

% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

msgAndWait(msgStr);

loopDiffMs = 1000*(loopNow-loopStart);
if loopDiffMs > msHold
    success = 1;
else
    success = 0;
end

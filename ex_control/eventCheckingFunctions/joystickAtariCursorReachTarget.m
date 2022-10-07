function success = joystickAtariCursorReachTarget(loopStart, loopNow, targX,targY,targRadius,startCursX, startCursY, cursorR, joystickDevInd, cursorColor, targWinCursRad, joystickPxPerSec)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes
persistent cursorPos loopOld
if isempty(loopOld)
    loopOld = loopStart;
end

if isempty(cursorPos)
    cursorPos = [startCursX, startCursY];
end

pixBoxLimit = 300;
if any(abs(cursorPos)>pixBoxLimit)
    cursorPos(startPos<-pixBoxLimit) = -(pixBoxLimit-1);
    cursorPos(startPos>pixBoxLimit) = pixBoxLimit-1;
end

% for the control comptuer to track where target is
yellow  = [255 255 0];
winColors = yellow;
drawFixationWindows(targX,targY,targWinCursRad,winColors);

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];


loopDiff = loopNow - loopOld;
loopOld = loopNow;
joystickMaxValue = 32768;
adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;

% compute the new cursor position (and don't let it get out of the bounding
% box)
posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
posMv = posMv.*[1 -1]; % invert y axis

cursorPos = cursorPos + posMv*adjFactorForVel;

if any(abs(cursorPos)>pixBoxLimit)
    % return to old position if it's out of the bounding box
    cursorPos = cursorPos - posMv*adjFactorForVel;
end

% record cursor position to MAT file and NEV file *if* the cursor has moved
if ~all(posMv==0)
    posShift = posShiftForCode - cursorPos;
    posShiftX = posShift(1);
    posShiftY = posShift(2);
    sendCode(codes.CURSOR_POS);
    sendCode(posShiftX);
    sendCode(posShiftY);
end

% compute how close the cursor is to the target
relPos = ([targX targY] - cursorPos);
distToTarget = sqrt(sum(relPos.^2));

switch size(targRadius,1)
    case 1 %circular window
        success = distToTarget < targWinCursRad;
    case 2 %rectangular window
        success = all(abs(distToTarget)<abs(targRadius),1);
    otherwise
        error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
end

if keyboardEvents()
    success = -1;
end

% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

msgAndWait(msgStr);

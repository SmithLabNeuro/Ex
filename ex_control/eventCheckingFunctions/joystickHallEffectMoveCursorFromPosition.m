function [success, msgStr, fixWinOutput] = joystickHallEffectMoveCursorFromPosition(~, ~, positionHold, distanceTolerance, pixelDistForMaxJoystickPos, cursorObjId, cursorR, cursorColor, immediateFailOnReturn)
% success if the cursor leaves positionHold
% failure if the cursor *doesn't* leave positionHold

global codes

pixBoxLimit = pixelDistForMaxJoystickPos;

if nargin < 9
    immediateFailOnReturn = false;
end

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% if focusFlag
%     % focus object must be drawn
%     focusColorDisp = [focusColor(1) focusColor(2) focusColor(3)];
% end

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

% compute the new cursor position (and don't let it get out of the bounding
% box)
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
xVal = xVal*pixBoxLimit;
yVal = yVal*pixBoxLimit;
cursorPos = [xVal, yVal]; 

posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);


% compute how close the cursor is to the target


% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('set %i oval 0 %i %i %i %i %i %i', [cursorObjId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
positionXHold = positionHold(1);
positionYHold = positionHold(2);
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xyLargeEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, false);

if xyLargeEnough
    success = 1;
else
    if immediateFailOnReturn
        success = -1;
    else
        success = 0;
    end
end

cursorR = 5;
yellow  = [255 255 0];
winColors = yellow;

fixWinOutput = {[positionXHold cursorPosDisp(1)], [positionYHold cursorPosDisp(2)], [distanceTolerance cursorR],winColors};

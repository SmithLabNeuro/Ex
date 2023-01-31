function [success, msgStr, fixWinOutput, extraOutput] = joystickHallEffectCursorReachOneOfMultipleTargets(~,~, targsX,targsY, targsRewardIndices, cursorObjectId, cursorR, cursorColor, targWinCursRad, pixelDistForMaxJoystickPos)
% success if the cursor reaches any of the targets
% failure if the cursor *doesn't* reach a target.
% failure if the cursor reaches the wrong target.

global codes

pixBoxLimit = pixelDistForMaxJoystickPos;

% for the control computer to track where target is
yellow  = [255 255 0];
winColors = yellow;

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
% I think this posShift, a double, gets cast to an int in unixSendByte,
% which sendCode calls... keep in mind in case weird things appear, but
% it's been working until now...
posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);

numTargs = length(targsX);
distToEachTarget = zeros(numTargs , 1);
% Compute distance of cursor to each target
for targetIdx = 1:numTargs
    currTargX = targsX{targetIdx};
    currTargY = targsY{targetIdx};
    relPos = ([currTargX currTargY] - cursorPos);
    distToEachTarget(targetIdx) = sqrt(sum(relPos.^2));
end

% Check if distance to any Target is less than targWinCursRad
success = any(distToEachTarget < targWinCursRad);
% if reach is successful, identify which target animal chose
if success
    % Identify index of target it is closest to for success (1)
    [~, minIdx] = max(distToEachTarget < targWinCursRad);
    % Report back the reward index that he selected
    extraOutput = minIdx;
else
    extraOutput = nan;
end

% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

% msgAndWait(msgStr);


% drawFixationWindows([targX cursorPosDisp(1)],[targY cursorPosDisp(2)],[targWinCursRad cursorR],winColors);
fixWinOutput = {[targsX{1} cursorPosDisp(1)], [targsY{1} cursorPosDisp(2)], [targWinCursRad cursorR], winColors};

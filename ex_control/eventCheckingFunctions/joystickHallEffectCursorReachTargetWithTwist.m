function [success, msgStr, fixWinOutput, extraOutput] = joystickHallEffectCursorReachTargetWithTwist(~,~, targX,targY,targRadius, cursorObjectId, cursorR, cursorColor, targWinCursRad, pixelDistForMaxJoystickPos, notTargX, notTargY, distTol, angleZHold, angleTolerance)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes

pixBoxLimit = pixelDistForMaxJoystickPos;

% for the control computer to track where target is
yellow  = [255 255 0];
winColors = yellow;
% drawFixationWindows(targX,targY,targWinCursRad,winColors);

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% these are the 0 values around which we'll d    [reachTarget, funcSuccesses] = waitForEvent(movementTime, {@joystickHallEffectCursorReachTarget, @joystickAtariHold, @turnOffTargetAfterMs}, {movementToTargetJoystickHEArgs, joystickAtariHoldArgs, turnOffTargetAfterMsArgs});
etermine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

% compute the new cursor position (and don't let it get out of the bounding
% box)
[xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();
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

% compute how close the cursor is to the target
relPos = ([targX targY] - cursorPos);
distToTarget = sqrt(sum(relPos.^2));

if ~isempty(notTargX)
    projOnBadTarg = [notTargX notTargY] * cursorPos';
    projOnGoodTarg = [targX targY] * cursorPos';
    distFromCent = sqrt(sum(cursorPos.^2));
    targetDist = sqrt(sum(targX^2 + targY^2));
    relPosFromBadTarg = [notTargX notTargY] - cursorPos;
    distToBadTarg = sqrt(sum(relPosFromBadTarg.^2, 2));
end

% see if joystick is still twisted:
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);
extraOutput = abs(zValAng);

switch size(targRadius,1)
    case 1 %circular window
        success = distToTarget < targWinCursRad;
        if ~isempty(notTargX) && ((any(projOnBadTarg>projOnGoodTarg) && (distFromCent > targetDist*distTol)) || (projOnGoodTarg<0 && distFromCent>targetDist*distTol)) %qany(distToBadTarg < targWinCursRad)
            success = -1;
        end
        if ~angLargeEnough
            % HE not twisted
            success = -1;
        end
%         if targY==0 && targX > 0 && distToTarget < targWinCursRad
%             disp(distToTarget)
%         end
    case 2 %rectangular window
        success = all(abs(distToTarget)<abs(targRadius),1);
        if ~isempty(notTargX) && any(all(abs(distToBadTarg)<abs(targRadius'),2))
            success = -1;
        end
    otherwise
        error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
end

% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

% msgAndWait(msgStr);


% drawFixationWindows([targX cursorPosDisp(1)],[targY cursorPosDisp(2)],[targWinCursRad cursorR],winColors);
fixWinOutput = {[targX cursorPosDisp(1)], [targY cursorPosDisp(2)], [targWinCursRad cursorR], winColors};

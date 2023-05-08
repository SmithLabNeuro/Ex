function [success, msgStr, fixWinOutput, funcOutputs] = joystickHallEffectHoldTargetForMs(loopStart, loopNow, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos, msHold, cursorObjectId, cursorR, cursorColorDisp)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position
global codes
persistent holdTargStart

if isempty(holdTargStart)
    holdTargStart = loopStart;
end

successUnlessFail = false;
success = 0;
msgStr = '';
[xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();
pixBoxLimit = pixelDistForMaxJoystickPos;
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% angLargeEnough = zValAng > (angleZHold-angleTolerance);
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);

% prStr = sprintf('angle %6.2f\n', zValAng);
% prStrB = prStr(1:end);
% prStrB(:) = sprintf('\b');
% fprintf(prStrB)
% fprintf(prStr)

% if the hold time has passed, success no matter what (so this function
% doesn't care what you do past the hold time)
% loopDiffMs = 1000*(loopNow-loopStart);
loopDiffMs = 1000*(loopNow-holdTargStart);
if loopDiffMs > msHold
    success = 1;
else
    if ~(xySmallEnough) %~(xySmallEnough && angLargeEnough) % ~(xySmallEnough && buttonPress) %
        holdTargStart = loopNow; % reset loop? yes!
%         success = -1;
    end
end

yellow  = [255 255 0];
winColors = yellow;
cursorPos = [xVal, yVal]; 

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);

cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
cursorIndicatorAngle = 20; % degree wedge
cursAngDispWedgeVals = [cursorR, cursorIndicatorAngle/2, zValAng]';

% draw cursor
msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

numWindows = 3;
maxSizeInfoVals = 3;
sizeInfo = nan(maxSizeInfoVals, numWindows);
sizeInfo(1:length(distanceTolerance),1) = distanceTolerance;
sizeInfo(1:length(cursorR),2) = cursorR;
sizeInfo(1:length(cursAngDispWedgeVals),3) = cursAngDispWedgeVals;
fixWinOutput = {[positionXHold cursorPosDisp(1) cursorPosDisp(1)], [positionYHold cursorPosDisp(2) cursorPosDisp(2)], sizeInfo,winColors};

funcOutputs = zValAng;
end

function passTolCheck = checkWithinTolerance(actual, expected, tolerance, checkInIfTrue)
    passTolCheck = false;
    % flag to see if the one must stay within the tolerance or be outside
    % the tolerance
    checkOutIfTrue = ~checkInIfTrue;
    if checkInIfTrue && (abs(actual - expected) < tolerance)
        passTolCheck = true;
    elseif checkOutIfTrue &&  (abs(actual-expected) > tolerance)
        passTolCheck = true;
    end
end



% function [success, msgStr, fixWinOutput, funcOutputs] = joystickHallEffectHoldTargetForMs(loopStart, loopNow, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos, msHold, cursorObjectId, cursorR, cursorColorDisp)
% % success if the *joystick* (not the cursor) reaches correct position
% % failure if the joystick does not reach correct position
% global codes
% persistent holdTargStart
% 
% if isempty(holdTargStart)
%     holdTargStart = loopStart;
% end
% 
% successUnlessFail = false;
% success = 0;
% msgStr = '';
% [xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();
% pixBoxLimit = pixelDistForMaxJoystickPos;
% xVal = xVal * pixBoxLimit;
% yVal = yVal * pixBoxLimit;
% 
% % check that x, y, and z position are at desired hold positions
% distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
% xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% % angLargeEnough = zValAng > (angleZHold-angleTolerance);
% angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);
% 
% % prStr = sprintf('angle %6.2f\n', zValAng);
% % prStrB = prStr(1:end);
% % prStrB(:) = sprintf('\b');
% % fprintf(prStrB)
% % fprintf(prStr)
% 
% % if the hold time has passed, success no matter what (so this function
% % doesn't care what you do past the hold time)
% % loopDiffMs = 1000*(loopNow-loopStart);
% loopDiffMs = 1000*(loopNow-holdTargStart);
% if loopDiffMs > msHold
%     success = 1;
% else
%     if ~(xySmallEnough) %~(xySmallEnough && angLargeEnough) % ~(xySmallEnough && buttonPress) %
%         holdTargStart = loopNow; % reset loop? yes!
% %         success = -1;
%     end
% end
% 
% yellow  = [255 255 0];
% winColors = yellow;
% cursorPos = [xVal, yVal]; 
% 
% % these are the 0 values around which we'll determine cursor position (for
% % saving to the NEV file and also the MAT file)
% posX0 = 10000;
% posY0 = 10000;
% posShiftForCode = [posX0 posY0];
% 
% posShift = posShiftForCode + cursorPos;
% posShiftX = posShift(1);
% posShiftY = posShift(2);
% sendCode(codes.CURSOR_POS);
% sendCode(posShiftX);
% sendCode(posShiftY);
% 
% cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
% cursorIndicatorAngle = 20; % degree wedge
% cursAngDispWedgeVals = [cursorR, cursorIndicatorAngle/2, zValAng]';
% 
% % draw cursor
% msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
% 
% numWindows = 3;
% maxSizeInfoVals = 3;
% sizeInfo = nan(maxSizeInfoVals, numWindows);
% sizeInfo(1:length(distanceTolerance),1) = distanceTolerance;
% sizeInfo(1:length(cursorR),2) = cursorR;
% sizeInfo(1:length(cursAngDispWedgeVals),3) = cursAngDispWedgeVals;
% fixWinOutput = {[positionXHold cursorPosDisp(1) cursorPosDisp(1)], [positionYHold cursorPosDisp(2) cursorPosDisp(2)], sizeInfo,winColors};
% 
% funcOutputs = zValAng;
% end
% 
% function passTolCheck = checkWithinTolerance(actual, expected, tolerance, checkInIfTrue)
%     passTolCheck = false;
%     % flag to see if the one must stay within the tolerance or be outside
%     % the tolerance
%     checkOutIfTrue = ~checkInIfTrue;
%     if checkInIfTrue && (abs(actual - expected) < tolerance)
%         passTolCheck = true;
%     elseif checkOutIfTrue &&  (abs(actual-expected) > tolerance)
%         passTolCheck = true;
%     end
% end
% 
% 

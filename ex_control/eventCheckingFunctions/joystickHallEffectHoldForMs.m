function success = joystickHallEffectHoldForMs(loopStart, loopNow, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, msHold)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position

success = 0;
[xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();
pixBoxLimit = 300;
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% angLargeEnough = zValAng > (angleZHold-angleTolerance);
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);

% check for a position change
if ~(xySmallEnough && angLargeEnough)
    success = -1;
else
    loopDiffMs = 1000*(loopNow-loopStart);
    if loopDiffMs > msHold
        success = 1;
    else
        success = 0;
    end
end

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
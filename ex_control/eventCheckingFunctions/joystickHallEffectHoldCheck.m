function success = joystickHallEffectHoldCheck(~, ~, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position

success = 0;
[xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% angLargeEnough = zValAng > (angleZHold-angleTolerance);
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);

if ~(xySmallEnough && angLargeEnough)
    success = -1;
else
    success = 1;
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
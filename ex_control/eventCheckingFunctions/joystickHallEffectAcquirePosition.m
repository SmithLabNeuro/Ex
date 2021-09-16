function success = joystickHallEffectAcquirePosition(~, ~, joystickButtonX, joystickButtonY, joystickDevInd, checkAll, joystickOnMsgs, joystickOffMsgs)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position
sendMessages = false; 
if nargin < 6
    % assumes checkAll is not provided and sets default value
    checkAll = true;
elseif nargin == 7
    error('Missing either joystick ON messages or joystick OFF messages.')
else
    sendMessages = true;
end

successCheck = zeros(length(joystickDevInd), 1);

[xVal, yVal, zValAng, button] = sampleHallEffectJoystick()
raise('Haven''t completed this code yet')
% check for position being correct
if ~sum(abs(joystickXY - [joystickButtonX, joystickButtonY]))
    successCheck(i) = 1;
    if sendMessages
        msgAndWait(joystickOnMsgs{i});
    end
else
    successCheck(i) = 0;
    if sendMessages
        msgAndWait(joystickOffMsgs{i});
    end
end

if checkAll
    success = all(successCheck);
else
    success = any(successCheck);
end

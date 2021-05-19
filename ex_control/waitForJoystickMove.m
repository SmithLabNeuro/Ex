function [success, cursorPos] = waitForJoystickMove(waitTime, cursorR, joystickDevInd, joystickPxPerSec, cursorColor, varargin)

pixBoxLimit = 300;
if ~isempty(varargin)
    startPos = round(varargin{1});
    if any(abs(startPos)>pixBoxLimit)
        startPos(startPos<-pixBoxLimit) = -(pixBoxLimit-1);
        startPos(startPos>pixBoxLimit) = pixBoxLimit-1;
    end
else
    startPos = [0, 0];
end


startX = 0;
startY = 0;
posCurr = [startX startY];
posCurrFirst = posCurr;
cursorPos = posCurrFirst;
thisStartFirst = tic;
loopTop = GetSecs; % just initialize for the first round...
joystickMaxValue = 32768;
while toc(thisStartFirst)*1000 <= waitTime
    loopTopPrev = loopTop;
    loopTop = GetSecs;
    
    
    loopDiff = loopTop - loopTopPrev; % this is loop time in *seconds*
    adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;
    
    posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
    posMv = posMv.*[1 -1]; % invert y axis
    posCurr = posCurr + posMv*adjFactorForVel;
    if ~isequal(posCurr, posCurrFirst)
        break
    end
end

if isequal(posCurr, posCurrFirst)
    success = false;
    return
else
    success = true;
    cursorPos = posCurr;
end
function trialSuccess = waitForMSJoystick(waitTime,fixX,fixY,r, joystickDevInd, cursorColor, cursorR, joystickPxPerSec)
% function success = waitForMSJoystick(waitTime,fixX,fixY,r)
% 
% ex trial helper function: waits for t ms, checking to ensure that the
% joystick cursor remains within the joystick target window.  If time
% expires, trialSuccess = 1, but if fixation is broken first, trialSuccess
% returns 0
%
% waitTime: time to maintain fixation (in ms)
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window

global params;

    winColors = [255 255 0];
   
        
    if nargin < 8
        error('waitForMSJoystick can have exactly 8 input arguments');
    end
    
    drawFixationWindows(fixX,fixY,r,winColors);
    % draw the cursor at the fixation window...
    msg('set 3 oval 0 %i %i %i %i %i %i', [fixX fixY cursorR cursorColor(1) cursorColor(2) cursorColor(3)]);
    posCurr = [fixX fixY];
    
    trialSuccess = 1;
    thisStart = tic;

    joystickMaxValue = 32768;
    loopTop = GetSecs;

    while (toc(thisStart)*1000) <= waitTime
        
        loopTopPrev = loopTop;
        loopTop = GetSecs;
                
        
        loopDiff = loopTop - loopTopPrev; % this is loop time in *seconds*
        adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;
        
        posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
        posCurr = posCurr + posMv*adjFactorForVel;
        %           eyePos = eyePos - [fixX fixY];
        
        relPos = ([fixX fixY] - posCurr);
        relPosPix = sqrt(sum(relPos.^2));
    
        switch size(r,1)
            case 1, %circular window
                inWin = sum(relPos.^2)<r.^2;
            case 2, %rectangular window
                inWin = all(abs(relPos)<abs(r),1);
            otherwise
                error('EX:waitForMS:badRadius','Radius must have exactly 1 or 2 rows');
        end;
        
        
        if keyboardEvents()|| ~inWin
            trialSuccess = 0;
            break;
        end
%         if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
        
    end
        
    if nargin > 2
        drawFixationWindows()
    end
end

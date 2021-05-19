function trialSuccess = waitForMSJoystickAndFix(waitTime,fixX,fixY,fixR,cursX, cursY, cursR, joystickDevInd, cursorColor, cursorR, joystickPxPerSec)
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
   
        
    if nargin < 9
        error('waitForMSJoystick can have exactly 8 input arguments');
    end
    
    drawFixationWindows(cursX,cursY,cursR,winColors);
    % draw the cursor at the fixation window...
    msg('set 3 oval 0 %i %i %i %i %i %i', [cursX cursY cursorR cursorColor(1) cursorColor(2) cursorColor(3)]);
    posCurr = [cursX cursY];
    
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
        
        relPos = ([cursX cursY] - posCurr);
        relPosPix = sqrt(sum(relPos.^2));
    
        switch size(cursR,1)
            case 1, %circular window
                inWinJoy = sum(relPos.^2)<cursR.^2;
            case 2, %rectangular window
                inWinJoy = all(abs(relPos)<abs(cursR),1);
            otherwise
                error('EX:waitForMSJoystickAndFix:badRadius','Radius must have exactly 1 or 2 rows');
        end;
        
        d=samp;
            eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013

 %           eyePos = eyePos - [fixX fixY];

            relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window        
            switch size(fixR,1)
                case 1, %circular window        
                    inWin = sum(relPos.^2,1)<fixR.^2; 
                case 2, %rectangular window
                    inWin = all(abs(relPos)<abs(fixR),1);
                otherwise
                    error('EX:waitForMSJoystickAndFix:badRadius','Radius must have exactly 1 or 2 rows');
            end;
        
        
        if keyboardEvents()|| ~inWin || ~inWinJoy
            trialSuccess = 0;
            break;
        end
%         if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
        
    end
        
    if nargin > 2
        drawFixationWindows()
    end
end

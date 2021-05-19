function [trialSuccess, reason, cursorPos] = waitForJoystickWhileFix(waitTime,fixX,fixY,fixR,targX, targY, targWinCursRad, joystickDevInd, cursorColor, cursorMoveColor, cursorR, joystickPxPerSec, varargin)
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
global codes


pixBoxLimit = 300;
if ~isempty(varargin);
    startPos = round(varargin{1});
    if any(abs(startPos)>pixBoxLimit)
        startPos(startPos<-pixBoxLimit) = -(pixBoxLimit-1);
        startPos(startPos>pixBoxLimit) = pixBoxLimit-1;
    end
else
    startPos = [0, 0];
end

cursorColorDisp = cursorColor;
breakFix = 0;
    winColors = [255 255 0];
    
% these are the 0 values around which we'll determine cursor position
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];


   
        
    if nargin < 9
        error('waitForMSJoystick can have exactly 8 input arguments');
    end
    
    startX = startPos(1);
    startY = startPos(2);

    drawFixationWindows([targX fixX],[targY fixY],[targWinCursRad fixR],winColors);
    % draw the cursor at the fixation window...

%     msg('set 3 oval 0 %i %i %i %i %i %i', [startX startY cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
%     msgAndWait('obj_on 3');
%     sendCode(codes.CURSOR_ON)
    cursOn = true;

    
        posCurr = [startX startY];
        posCurrFirst = posCurr;
        cursorPos = posCurrFirst;
    
    trialSuccess = 0;
    joystickMaxValue = 32768;
    
    thisStartFirst = tic;
    
    loopTop = GetSecs;
    
%     while toc(thisStartFirst)*1000 <= 500
%         loopTopPrev = loopTop;
%         loopTop = GetSecs;
%                 
%         
%         loopDiff = loopTop - loopTopPrev; % this is loop time in *seconds*
%         adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;
%         
%         posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
%         posMv = posMv.*[1 -1]; % invert y axis
%         posCurr = posCurr + posMv*adjFactorForVel;
%         if ~isequal(posCurr, posCurrFirst)
%             break
%         end
%     end
%     
%     if isequal(posCurr, posCurrFirst)
%         reason = 'IGNORED';
%         return
%     end
    
    thisStart = tic;

    loopTop = GetSecs;
    choice = 0;

    while (toc(thisStart)*1000) <= waitTime && ~choice
        
        loopTopPrev = loopTop;
        loopTop = GetSecs;
                
        
        loopDiff = loopTop - loopTopPrev; % this is loop time in *seconds*
        adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;
        
        posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
        posMv = posMv.*[1 -1]; % invert y axis
        posPrev = posCurr;
        posCurr = posCurr + posMv*adjFactorForVel;
        %           eyePos = eyePos - [fixX fixY];
        
        relPos = ([targX targY] - posCurr);
        relPosPix = sqrt(sum(relPos.^2));
        cursorPos =  posCurr;

        if any(abs(cursorPos)>pixBoxLimit)
            posNew = posOld;
            cursorPos = posNew; %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013 %-no longer supports polynomials with order>1, -acs09dec2015
        end

        if sum(sqrt(cursorPos.^2))>0%mod(floor(sum(sqrt(cursorPos.^2))/50), 2) %~mod(floor(toc(thisStart)*5)+1, 2)
            if ~isequal(cursorMoveColor, cursorColorDisp)
                sendCode(codes.CURSOR_OFF)
                cursOn = false;
            end
            cursorColorDisp = cursorMoveColor;
        else
            if ~isequal(cursorColor, cursorColorDisp)
                sendCode(codes.CURSOR_ON)
                cursOn = true;
            end
            cursorColorDisp = cursorColor;
        end
        
        cursorPosDisp = round(cursorPos); % rounding corrects %i issues... (an integer must be passed!)
        msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
        
        msgAndWait(msgStr);
        
        if ~all(posPrev == posCurr)
            
            %         sccom = [codes.CURSOR_POS posNew-posShiftForCode]'
            % recognize there's an inversion, but posShift should be
            % *positive* to help it get sent as an unsigned INT in sendCode
            % (otherwise it sends (value + UINT_MAX+1) to the NEV
            posShift = posShiftForCode - posCurr; 
            posShiftX = posShift(1);
            posShiftY = posShift(2);
            sendCode(codes.CURSOR_POS);
            sendCode(posShiftX);
            sendCode(posShiftY);
        end
    
        switch size(targWinCursRad,1)
            case 1, %circular window
                inWinJoy = sum(relPos.^2)<targWinCursRad.^2;
            case 2, %rectangular window
                inWinJoy = all(abs(relPos)<abs(targWinCursRad),1);
            otherwise
                error('EX:waitForMSJoystickAndFix:badRadius','Radius must have exactly 1 or 2 rows');
        end;
            choice = find([true,inWinJoy],1,'last')-1; %note in the case that the gaze is inside two or more overlapping windows, the choice is the last window provided. -acs09dec2015
            if choice
                
                trialSuccess = 1;
                break;
            end
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
        
        
        if keyboardEvents()|| ~inWin
            trialSuccess = 0;
            breakFix = 1;
            break;
        end
%         if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
%             drawFixationWindows([targX ,[targY fixY],[targWinCursRad fixR]);%,winColors);
drawFixationWindows([targX cursorPosDisp(1) fixX],[targY cursorPosDisp(2) fixY],[targWinCursRad cursorR fixR] ,winColors);

    end
    
     cursorColorDisp = cursorColor;
     if ~breakFix
%         if ~choice
%             pause(0.5)
        if ~isequal(cursorMoveColor, cursorColor)
            msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
            msgAndWait(msgStr);
            if ~cursOn
                sendCode(codes.CURSOR_ON)
            end
            pause(0.1)
        end
     end
    
     if ~trialSuccess
         if breakFix
             reason = 'BROKE_FIX';
         else
             reason = 'WRONG_TARG';
         end
     else
         reason = 'CORRECT';
     end
     
        
    if nargin > 2
        drawFixationWindows()
    end
    
end

function [choice, reason, cursorPos] = waitForJoystick(waitTime,targX,targY,r, cursorR, joystickDevInd, cursorColor, cursorMoveColor, targWinCursRad, joystickPxPerSec, varargin)
% function choice = waitForFixation(waitTime,fixX,fixY,r)
% function choice = waitForFixation(waitTime,fixX,fixY,r,winColors)
%
% ex trial helper function: looks for eye positions within a window,
% returning either when the eye enters the window or when time expires
%
% waitTime: time before function failure (in ms)
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window
% winColors: If specified, provides a list of colors (N x 3) for the
% fixation windows drawn on the user screen
% stability: a threshold, in standard deviation units, of the last 5
% position samples. Below this limit the eye position is considered stable.
%
% returns a value "choice" that indicates the window that was fixated, or
% 0 if no window was fixated.

% revised 09dec2015 by Adam Snyder - removed call to KeyboardEvents inside
% the while loop, which was slowing things down. This means you can't quit
% the trial or give juice during a call to waitForFixation (the keystrokes
% are buffered and will be dealt with when waitForFixation returns).

% Revised 24Feb2012 by Adam Snyder (adam@adamcsnyder.com) --based on an
% existing Ex function.
%
% Revised 23Oct2012 -ACS

global codes


assert(nargin>=4,'waitForFixationChoice must have fixation windows specified');
numWindows = unique([length(targX) length(targY) size(r,2)]);
assert(numel(numWindows)==1,'Fixation window parameters X, Y and R must be same size');
yellow  = [255 255 0];

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
% if nargin > 4
%     winColors = varargin{1};
%     if isempty(winColors),
%         winColors = yellow;
%     end;
% else
    winColors = yellow;
% end;

drawFixationWindows(targX,targY,targWinCursRad,winColors);

% red = [255 0 0];
posNew = startPos;% [0,0];
posFix = [0,0];
fixRad = 5;
runLoop = true;
targOn = true;

posMvPrev = [0,0];
giveJuiceBump = true;
cursOn = true;

cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];
% if false
%     msg('set 3 oval 0 %i %i %i %i %i %i', [posNew(1) posNew(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
%     % msg('mset 4 oval 0 %i %i %i %i %i %i',[targX targY r cursorColor(1) cursorColor(2) cursorColor(3)]);
%     msgAndWait('obj_on 3');
% end

% while toc(pauseStart) <= secsPauseCursorOn
%     posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
%     posMv = posMv.*[1 -1]; % invert y axis
% 
%     
%     if ~isempty(posMvPrev) && ~isequal(posMv, posMvPrev)
% %         runLoop = false;
%         pause(0.5)
%         giveJuiceBump = false;
%         runLoop = false;
%         disp('BLAHBLAH')
%         disp(posMv)
%         break;
%     else
%         posMvPrev = posMv;
%     end
% end



thisStart = tic;
loopTop = GetSecs; % just initialize for the first round...

% these are the 0 values around which we'll determine cursor position
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

choice = 0;
% posMvPrev = [];
giveFirstJuice = true;
timeGotMoveJuice = 0;
joystickMaxValue = 32768;
loopDiffPrev = 0;
while (toc(thisStart)*1000)<=waitTime && choice<1 %changed from while 1 so that there are less commands to evaluate inside the loop
    
    loopTopPrev = loopTop;
    loopTop = GetSecs;
    loopDiff = loopTop - loopTopPrev; % this is loop time in *seconds*
    
%     if toc(thisStart) <= secsCursorNoMove
%         posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
%         posMv = posMv.*[1 -1]; % invert y axis
%         
%         if ~isempty(posMvPrev) && ~isequal(posMv, posMvPrev)
%             runLoop = false;
%             pause(0.5)
%             giveJuiceBump = false;
%             runLoop = false;
%             disp('BLAHBLAH')
%             disp(posMv)
%             break;
%         else
%             posMvPrev = posMv;
%         end
%     elseif giveJuiceBump
%                 msg('set 1 oval 0 %i %i %i %i %i %i',[posFix(1) posFix(2) fixRad 0 255 0]);
%                 % UNCOMMENT if you want MGR, but with no delay
%                 if targOn
%                     msgAndWait('obj_off 4');
%                     targOn = false;
%                 end
% %         giveJuice(1,1,200);
% %         giveJuiceBump = false;
% %         if loopDiffPrev
% %             loopDiff = loopDiffPrev;
% %         end
%     end
    
    
    adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;
    
    posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
    posMv = posMv.*[1 -1]; % invert y axis

    
    if ~isempty(posMvPrev) && ~isequal(posMv, posMvPrev)

            % UNCOMMENT for juice ever two joystick clicks
%         if mod(timeGotMoveJuice, 2)==0
%             % unfortunately unclear how to record *amount* of reward
%             choice = 1;
%             sendCode(codes.REWARD);
%             giveJuice(1,1,50);
%         end
%         timeGotMoveJuice = timeGotMoveJuice+1;
        
        giveClickJuice = true;
    else
        giveClickJuice = false;
    end
    
    posMvPrev = posMv;
    posOld = posNew;
    posNew = posNew + posMv*adjFactorForVel;
        
    cursorPos = posNew;

    if any(abs(cursorPos)>pixBoxLimit)
        posNew = posOld;
        cursorPos = posNew; %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013 %-no longer supports polynomials with order>1, -acs09dec2015
    end
    
    if ~all(posOld == posNew)
        
%     posMv
%     posNew
%         sccom = [codes.CURSOR_POS posNew-posShiftForCode]'
        posShift = posShiftForCode - posNew;
        posShiftX = posShift(1);
        posShiftY = posShift(2);
        sendCode(codes.CURSOR_POS);
        sendCode(posShiftX);
        sendCode(posShiftY);
    end

    
    
    
    relPos = ([targX targY] - cursorPos);
    relPosPix = sqrt(sum(relPos.^2));
    targDist = sqrt(sum([targX targY].^2));
%     if relPosPix < 0.5*targDist && giveFirstJuice
%         giveFirstJuice = false;
%         giveJuice(1,1,75);
%     end
    
    
%     if giveClickJuice
%         % adding + r gives a baseline juice around the target that
        % encompasses the initial center point; 
%         % the max ensures only positive TTLs
% %         relPosPix
% %         dv=-0.5*(relPosPix)/targDist
% %         tot=bestNonrewardTTL*(-0.5*(relPosPix)/targDist)+bestNonrewardTTL
%         juiceTTL = max(bestNonrewardTTL*(-0.5*(relPosPix)/targDist)+bestNonrewardTTL, 0);
%         
%         giveJuice(1,1,juiceTTL);
%         % unfortunately unclear how to record *amount* of reward
%         sendCode(codes.REWARD);
%     end
    
%     
%     disp(sum(relPos.^2))
%     disp(r.^2)
% if relPosPix <  targWinRad + cursorR*5.5 && targOn
%     msgAndWait('obj_off 4');
%                     targOn = false;
% end
    switch size(r,1)
        case 1 %circular window        
%                         inWin = relPosPix < targWinRad + cursorR*2.5;
            inWin = relPosPix < targWinCursRad;
        case 2 %rectangular window
            inWin = all(abs(relPos)<abs(r),1);
        otherwise
            error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
    end
    choice = find([true,inWin],1,'last')-1; %note in the case that the gaze is inside two or more overlapping windows, the choice is the last window provided. -acs09dec2015

    if keyboardEvents()
        choice = 0;
        break;
    end
    
    cursorPosDisp = round(cursorPos); % rounding corrects %i issues... (an integer must be passed!)
if sum(sqrt(cursorPos.^2))>0%mod(floor(sum(sqrt(cursorPos.^2))/50), 2) %~mod(floor(toc(thisStart)*5)+1, 2) 
    if ~isequal(cursorMoveColor, cursorColorDisp)
        sendCode(codes.CURSOR_OFF)
        cursOn = false;
    end
    cursorColorDisp = cursorMoveColor;
    msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
    
else
    if ~isequal(cursorMoveColor, cursorColorDisp)
        sendCode(codes.CURSOR_ON)
        cursOn = true;
    end
    msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
end
msgAndWait(msgStr);

drawFixationWindows([targX cursorPosDisp(1)],[targY cursorPosDisp(2)],[targWinCursRad cursorR],winColors);


    %     msgAndWait('mset 3 oval 0 %i %i %i %i %i %i', [0 0 cursorR cursorColor(1) cursorColor(2) cursorColor(3)]);
%     if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForFixation exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
    loopDiffPrev = loopDiff;
end

if ~choice
    reason = 'WRONG_TARG';
%     cursorColorDisp = cursorColor;
%     msgStr = sprintf('mset 3 oval 0 %i %i %i %i %i %i', [cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
%     msgAndWait(msgStr);
%     if ~choice
%         pause(0.5)
%     else
%         pause(0.5)
%     end
else
    reason = 'CORRECT';
end

% if choice == 1
%     timeJuice = 1;
%     giveJuice(timeJuice);
% end
% drawFixationWindows()

end


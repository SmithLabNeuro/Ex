function trialSuccess = waitForMS(waitTime,fixX,fixY,r,varargin)
% function success = waitForMS(waitTime,fixX,fixY,r)
% 
% ex trial helper function: waits for t ms, checking to ensure that the eye
% remains within the fixation window.  If time expires, trialSuccess = 1, 
% but if fixation is broken first, trialSuccess returns 0
%
% waitTime: time to maintain fixation (in ms)
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window
%
% 2015/08/14 by Adam Snyder and Matt Smith. Now aloopTopllows user to pass a
% 'recenterFlag' such that the fixX and fixY will be ignored and instead
% the current eye position is used.
%

global params;

    winColors = [255 255 0];
    recenterFlag = false;
    if nargin > 4        
        vx = 1;
        while vx <= numel(varargin),
            switch class(varargin{vx})
                case 'char'
                    recenterFlag = varargin{vx+1};      
                    vx = vx+2;
                otherwise
                    if ~isempty(varargin{vx}),
                        winColors = varargin{vx};
                    end;
                    vx = vx+1;
            end;
        end;   
    end;
    
    if recenterFlag,
        d=samp;
        eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013
        fixX = eyePos(1);
        fixY = eyePos(2); 
%        fixX = eyePos(1)+fixX; % removed the fixX/Y addition here - that's
%        not right % MAS Sept2019
%        fixY = eyePos(2)+fixY; %is the sign here correct? -ACS 14Aug2015
    end;
    
    if nargin >= 4
        drawFixationWindows(fixX,fixY,r,winColors);
    elseif nargin ~=1
        error('waitForMS can have exactly 1, 4 or 5 input arguments');
    end
    
    trialSuccess = 1;
    thisStart = tic;

    if nargin >= 4         
        while (toc(thisStart)*1000) <= waitTime
            loopTop = GetSecs;
            d=samp;
            eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013

 %           eyePos = eyePos - [fixX fixY];

            relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window        
            switch size(r,1)
                case 1, %circular window        
                    inWin = sum(relPos.^2,1)<r.^2; 
                case 2, %rectangular window
                    inWin = all(abs(relPos)<abs(r),1);
                otherwise
                    error('EX:waitForMS:badRadius','Radius must have exactly 1 or 2 rows');
            end;

            if keyboardEvents()||~inWin
                trialSuccess = 0;
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
            
        end
    else %don't worry about fixation window - this is essentially just a pause (can be broken with a key press)
        while (toc(thisStart)*1000) <= waitTime
            loopTop = GetSecs;            
            if keyboardEvents()
                trialSuccess = 0;
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
        end
    end
    if nargin > 2
        drawFixationWindows()
    end
end

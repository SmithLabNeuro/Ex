function choice = waitForFixation(waitTime,fixX,fixY,r,varargin)
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

global params;

assert(nargin>=4,'waitForFixationChoice must have fixation windows specified');
numWindows = unique([length(fixX) length(fixY) size(r,2)]);
assert(numel(numWindows)==1,'Fixation window parameters X, Y and R must be same size');
yellow  = [255 255 0];
if nargin > 4
    winColors = varargin{1};
    if isempty(winColors),
        winColors = yellow;
    end;
else
    winColors = yellow;
end;

drawFixationWindows(fixX,fixY,r,winColors);

thisStart = tic;

choice = 0;
while (toc(thisStart)*1000)<=waitTime && choice<1, %changed from while 1 so that there are less commands to evaluate inside the loop
    loopTop = GetSecs;
    d=samp;
    eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013 %-no longer supports polynomials with order>1, -acs09dec2015
    relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window        
    switch size(r,1)
        case 1, %circular window        
            inWin = sum(relPos.^2,1)<r.^2; 
        case 2, %rectangular window
            inWin = all(abs(relPos)<abs(r),1);
        otherwise
            error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
    end;
    choice = find([true,inWin],1,'last')-1; %note in the case that the gaze is inside two or more overlapping windows, the choice is the last window provided. -acs09dec2015

    if keyboardEvents()
        choice = 0;
        break;
    end
    if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForFixation exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
end
drawFixationWindows()

end


function [fixSuccess, freqSuccess, freqvals ]= waitForMS(waitTime,fixX,fixY,r,targetfreq, freqtolerance,tolerancefunction, targettime,juicedur,varargin)
% function success = waitForMS(waitTime,fixX,fixY,r)
% 
% ex trial helper function: waits for t ms, checking to ensure that the eye
% remains within the fixation window.  If time expires, fixSuccess = 1, 
% but if fixation is broken first, fixSuccess returns 0
%
% waitTime: time to maintain fixation (in ms)
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window
%
% 2015/08/14 by Adam Snyder and Matt Smith. Now allows user to pass a
% 'recenterFlag' such that the fixX and fixY will be ignored and instead
% the current eye position is used.
%
% 2016/5/24 Ryan Williamson added functionality to check on success based
% on digital codes

global params;



 % check inputs
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
        fixY = eyePos(2); % removed fixX/Y - MAS Sept 2019
    end;
    
    if nargin >= 4
        drawFixationWindows(fixX,fixY,r,winColors);
    elseif nargin ~=1
        error('waitForMS can have exactly 1, 4 or 5 input arguments');
    end
    
    fixSuccess = 1;
    
    thisStart = tic;
    
    if nargin >= 4 
      
         freqSuccess = 0;
         successclockflag = 0;
         freqrewardflag = 0;
         freqvals = [];
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
            % retrieve digital events from Ripple and check freq
            [count, timestamps, events] = xippmex('digin');
            events = double([events.parallel]');
            if count>0 
                freqeventinds = events>1000 & events<32000;
                if length(freqeventinds)>0
                    freqevents = (events(freqeventinds))-1000;
                    % save freqs in behavioral data
                    freqvals = [freqvals; freqevents];
                    
                    % check if freq is within threshold of target. For now this
                    % does not require the target to be held 
                    if ~isempty(find(tolerancefunction(freqevents,targetfreq,freqtolerance),1))%~isempty(find(abs(freqevents-targetfreq)<freqtolerance,1))
                        if successclockflag == 0
                            successclock = toc(thisStart)*1000;
                            successclockflag = 1;
                        end
                        if successclockflag && abs(successclock- toc(thisStart)*1000)>targettime
                      
                            abs(successclock- toc(thisStart)*1000)
                            freqSuccess = 1;
                           
                            successclockflag = 0;
                            giveJuice(1,juicedur,juicedur);
                        end
                    else
                        successclockflag = 0;
                    end
                end
            end
           
            if keyboardEvents()||~inWin
                fixSuccess = 0;
%                 if isempty(freqvals)
%                         freqvals = {freqvalstemp};
%                         here = 3
%                 else
%                         freqvals(end+1)={freqvalstemp};
%                         here = 4
%                 end
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
            
        end
%         if isempty(freqvals)
%             freqvals = {freqvalstemp};
%             here = 1
%         else
%             freqvals(end+1)={freqvalstemp};
%             here = 2
%         end
    else %don't worry about fixation window - this is essentially just a pause (can be broken with a key press)
        while (toc(thisStart)*1000) <= waitTime
            loopTop = GetSecs;            
            if keyboardEvents()
                fixSuccess = 0;
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
        end %waitTime
    end
    
    if nargin > 2
        drawFixationWindows()
    end              
end

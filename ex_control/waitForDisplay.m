function trialSuccess = waitForDisplay(fixX,fixY,r,objID)
% function trialSuccess = waitForDisplay(fixX,fixY,r,objID)
%
% ex trial helper function: waits for a message from the display
% computer, checking to ensure that the eye remains within the
% fixation window.  If the display responds, trialSuccess = 1, but if
% fixation is broken first, trialSuccess returns 0
%
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window
%
% ACS, 15mar2016: added argument objID, if supplied waits for that specific
% object to finish.

global params debug sockets;

drawFixationWindows(fixX,fixY,r);

trialSuccess = 1;
while true
    if matlabUDP2('check',sockets(1))
        [s1, s] = strtok(matlabUDP2('receive',sockets(1)));
        switch s1
            case 'done' %message received, don't need to wait any longer. -acs15mar2016
                if nargin<4||str2double(s)==objID,
                    break;
                else
                    if debug,
                        fprintf('waitForDisplay: object %s ended, but was waiting for object %s.\n', s, objID);
                    end;
                end;
            otherwise
                if debug,
                    fprintf('waitForDisplay: received message ''%s %s''.\n', s, s1);
                end;
        end;
    end;
    loopTop = GetSecs;
    d=samp;
    eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013

    relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window        
    switch size(r,1)
        case 1, %circular window        
            inWin = sum(relPos.^2,1)<r.^2; 
        case 2, %rectangular window
            inWin = all(abs(relPos)<abs(r),1);
        otherwise
            error('EX:waitForDisplay:badRadius','Radius must have exactly 1 or 2 rows');
    end;
    
    if keyboardEvents()||~inWin
        trialSuccess = 0;
        break;
    end
    if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForDisplay exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
end

drawFixationWindows();
end

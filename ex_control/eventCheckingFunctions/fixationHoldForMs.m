function [success, msgStr, fixWinOutput] = fixationHoldForMs(loopStart, loopNow, fixX, fixY, fixWinVals, msHold, successUnlessFail)

msgStr = '';
success = 0;
fixWinOutput = {};
purple  = [255 0 255];
winColors = purple;

if nargin<7
    successUnlessFail = false;
end

% if the hold time has passed, success no matter what
loopDiffMs = 1000*(loopNow-loopStart);
if loopDiffMs > msHold
    success = 1;
    
else
    eyeVoltage=samp;
    eyePos = projectCalibration(eyeVoltage(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013 %-no longer supports polynomials with order>1, -acs09dec2015
    relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window
    switch size(fixWinVals,1)
        case 1 %circular window
            eyesInWindow = sum(relPos.^2,1)<fixWinVals.^2;
            fixWinOutput = {fixX, fixY, fixWinVals, winColors};
        case 2 %rectangular window
            eyesInWindow = all(abs(relPos)<abs(fixWinVals),1);
            fixWinOutput = {};
        otherwise
            error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
    end
    
    if ~eyesInWindow
        % Don't penalize for samples being dropped
        if any(isnan(relPos))
            success = 0;
        else
            % eyes left the window - fail!
            success = -1;
        end
    else
        % success unless fail!
        if successUnlessFail
            success = 1;
        else
            % no success until time passed!
            success =  0;
        end
    end
end
function [success, msgStr, fixWinOutput] = fixationHold(~, ~, fixX, fixY, fixWinVals)

msgStr = '';
purple  = [255 0 255];
winColors = purple;

eyeVoltage=samp;
eyePos = projectCalibration(eyeVoltage(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013 %-no longer supports polynomials with order>1, -acs09dec2015
relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window
switch size(fixWinVals,1)
    case 1 %circular window
        success = sum(relPos.^2,1)<fixWinVals.^2;
        fixWinOutput = {fixX, fixY, fixWinVals, winColors};
    case 2 %rectangular window
        success = all(abs(relPos)<abs(fixWinVals),1);
        fixWinOutput = {};
    otherwise
        error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
end


if ~success
    % Don't penalize for samples being dropped
    if any(isnan(relPos))
        success = 1;
    else
        success = -1;
    
    end
end
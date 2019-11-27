function drawCalibration(n)
% function drawCalibration(n)
%
% used by the calibration procedure.  Draws up to n dots on the control
% screen.

global wins calibration;

if nargin<1, n = size(calibration{2},1); end; %if no input is passed, just draw all the dots. This helps eliminate magic numbers in runex. -ACS 03Sep2013
            
    for i = 1:n
        pt = calibration{2}(i,1:2); %changed to just use the first two points, because if polynomial regressors are used there might be more than two columns -ACS 29Oct2013
        pt = pt.*wins.pixelsPerVolt+wins.midV;
        Screen('FillOval',wins.voltageBG,wins.controlCalibDotColor,[pt - wins.controlCalibDotRad pt + wins.controlCalibDotRad]);
    end
end

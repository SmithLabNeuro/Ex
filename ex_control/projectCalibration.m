function calibrated = projectCalibration(in)
%Function to apply the calibration to the raw voltage values. -ACS
%29Oct2013

%removed nonlinear (polynomial and spherical) capabilities -acs 04dec2015

    global calibration

    %if params.getEyes,
    %    pos = pos .* wins.pixelsPerMV + wins.midV;
    %    pos(2) = wins.voltageDim(4) - pos(2); % flip Y coordinate
    %end;
    
    calibrated(:,1) = [in ones(size(in,1),1)] * calibration{3};
    calibrated(:,2) = [in ones(size(in,1),1)] * calibration{4};

end

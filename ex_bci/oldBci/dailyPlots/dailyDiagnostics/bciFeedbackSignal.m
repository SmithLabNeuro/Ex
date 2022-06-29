function [ fSig, popVecMag, popVecAngle ] = bciFeedbackSignal( weights, angles, trueAngle )
% weights - weights of each angle for each sample (nSamples x nAngles)
% angles  - angle in degrees (nAngles x 1)

    radAngles = deg2rad(angles);
    radTrueAngle = deg2rad(trueAngle);
    
    % compute the weighted vector and it's amplitude
    popVecAngle = weights*exp(radAngles.*1j);
    % compute angular difference with the true angle
    fSig = abs(angle( popVecAngle ./ exp(radTrueAngle.*1j) ));
    fSig = rad2deg(fSig);
    % compute population vector magnitude
    popVecMag = abs(popVecAngle);
    
end


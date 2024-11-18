function newVelocity = kalmanVelocityDemoDecoder(meanSpikeCount, currVelocity, modelParams, expParams)

M0 = modelParams.M0;
M1 = modelParams.M1;
M2 = modelParams.M2;
if isfield(modelParams, 'rotMat')
    rotMat = modelParams.rotMat;
else
    rotMat = eye(2);
end
newVelocity = rotMat*(M0 + M1*rotMat'*currVelocity + M2 * meanSpikeCount);
% xoutBase = repmat('0', 8, 8);
% bl = dec2bin(typecast(meanSpikeCount(1), 'uint8'));
% xoutBase(:, 8-size(bl, 2)+1:8) = bl;
% xout = double(xoutBase(:))-48;
% for i = 2:length(meanSpikeCount)
%     xoutBase = repmat('0', 8, 8);
%     bl = dec2bin(typecast(meanSpikeCount(i), 'uint8'));
%     xoutBase(:, 8-size(bl, 2)+1:8) = bl;
%     xout = xor(double(xoutBase(:))-48, xout);
% end
% find(xout)'
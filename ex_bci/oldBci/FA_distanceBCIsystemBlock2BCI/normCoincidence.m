function [ C ] = normCoincidence( X )
% Input:
%    - X: spike train matrix (nChannels x nT)
%
% Output:
%    - C: coincidences normalized by geometric mean
%
    
    s = X*X';
    normMat = diag(1./(sqrt(diag(s))));
    C = normMat*s*normMat;
    
    C(tril(true(size(C)))) = nan;
    
end


function [r, lagVec] = myXcorr(X, maxLag, verbose)
%MYXCORR Summary of this function goes here
%   Detailed explanation goes here

    if nargin<3
        verbose = false;
    end
    nNeurons = size(X,1);
    
    r = nan(nNeurons,nNeurons,2*maxLag+1);
    
    for ii=1:nNeurons
        if verbose
            fprintf('      Getting ccgs for neuron %d of %d...\n',ii,nNeurons);
        end
        for jj=(ii+1):nNeurons
            [tmp,lagVec] = xcorr(X(ii,:),X(jj,:),maxLag,'coeff');
            r(ii,jj,:) = tmp;
        end
    end
    
end


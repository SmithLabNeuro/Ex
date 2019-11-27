function out = deltaBCIDistMetric(latents1,latents2,mu1,mu2)
% distance metric for delta 2-bci experiment
    distance = @(x,mu)((sum((x-repmat(mu,1,size(x,2))).^2,1)));
    distance1 = distance(latents1,mu1);
    distance2 = distance(latents2,mu2);
    out = distance1-distance2;
end

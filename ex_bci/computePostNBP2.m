function posterior = computePostNBP2(counts,means)
    loglike = sum(-means + repmat(counts,1,size(means,2)).*log(means),1);
    test =exp(loglike-max(loglike));
    posterior=test./repmat(sum(test,2),1,length(test));
end
function posterior = computePostNBP(counts,means)
    logP = sum(-means + repmat(counts,1,size(means,2)).*log(means),1);
    const = max(logP);
    tmp = logP - const;
    tmp2 = log(sum(exp(tmp)))+const;
    logPost = logP - tmp2;
    posterior = exp(logPost);
    %posterior=test./repmat(sum(test,2),1,length(test));
end
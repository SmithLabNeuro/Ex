function [newsmoothlatent] = updateLatent(curcount,prevsmoothlatent,modelparams)

% numsteps = number of newton updates to run (10 should be good)
%  Ryan Williamson, 10/11/2018

    estParams = modelparams.estParams;
    z = fastfa_estep(curcount,estParams);
    if isempty(prevsmoothlatent)
        prevsmoothlatent = z.mean;
    end
    newsmoothlatent = modelparams.alpha*z.mean+(1-modelparams.alpha)*prevsmoothlatent;
end
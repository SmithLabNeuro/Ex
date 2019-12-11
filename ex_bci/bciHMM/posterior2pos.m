function [x,y] = posterior2pos(targ_angs,postProbs)
%
% targ_angs: all possible target angles (nStates x 1)
% postProbs: posteriors for each corresponding angle (nStates x 1)

    rad_vals = deg2rad(targ_angs);
    comp_vals = exp(rad_vals*1i);
    
    pos_comp = comp_vals'*postProbs;
    
    x = real(pos_comp);
    y = imag(pos_comp);
    
end


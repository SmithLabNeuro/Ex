function [pa,pdA,pdB] = prinangle(A,B)
    if nargout == 1
        sv = svd(orth(A)' * orth(B));
        pa = acos(sv);
    else
        Aorth = orth(A);
        Borth = orth(B);
        
        [UU, DD, VV] = svd(Aorth' * Borth, 'econ');
        
        pa = acos(diag(DD));
        pdA = Aorth * UU;
        pdB = Borth * VV;
    end
end
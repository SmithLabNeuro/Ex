function z_new = kf_filter(x,kf_params,z)

    d = kf_params.d;
    
    pival = kf_params.pival;
    V = kf_params.V;
    
    A = kf_params.A;
    Q = kf_params.Q;
    
    C = kf_params.C;
    R = kf_params.R;

    if isempty(z)
        mu_1_0 = pival;
        sig_1_0 = V;
    else
        mu_1_1 = z.mu;
        sig_1_1 = z.sig;
        mu_1_0 = A*mu_1_1;
        sig_1_0 = A*sig_1_1*A' + Q;
    end

    xt = x - d;

    kt = sig_1_0*C'*inv(C*sig_1_0*C'+R);
    z_new.mu = mu_1_0 + kt*(xt - C*mu_1_0);
    z_new.sig = sig_1_0 - kt*C*sig_1_0;

end


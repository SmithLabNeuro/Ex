function z_new = kf_onlineFilter(x,z,modelparams)
% x - spike counts
% z - previous latent
% modelparams - bci decoding parameters
    
    % project into FA space
    lowd_x = fastfa_estep(x,modelparams.fa_params);
    lowd_x = orthogonalize(lowd_x.mean,modelparams.fa_params.L);

    % apply KF
    kf_params = modelparams.kf_params;
    
    d = kf_params.d;
    
    A = kf_params.A;
    Q = kf_params.Q;
    C = kf_params.C;
    R = kf_params.R;

    if isempty(z)
        mu_1_0 = kf_params.pival;
        sig_1_0 = kf_params.V;
    else
        mu_1_1 = z.mu;
        sig_1_1 = z.sig;
        mu_1_0 = A*mu_1_1;
        sig_1_0 = A*sig_1_1*A' + Q;
    end

    xt = lowd_x - d;

    kt = sig_1_0*C'*inv(C*sig_1_0*C'+R);
    z_new.mu = mu_1_0 + kt*(xt - C*mu_1_0);
    z_new.sig = sig_1_0 - kt*C*sig_1_0;
    
end


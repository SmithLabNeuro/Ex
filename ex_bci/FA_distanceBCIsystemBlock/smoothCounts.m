function [meanvarnew] = smoothCounts(observedval,KFparams,meanvar)
meanvarnew.mu = 0;
meanvarnew.sig = 0;
%% Initialize parameter values
mutmin = meanvar.mu;
sigtmin = meanvar.sig;
A  = KFparams.A;
Q = KFparams.Q;
C = KFparams.C(:,1:end-1);
d = KFparams.C(:,end);
R = KFparams.R;

%% One-step prediction
muonestep = A*mutmin;
sigonestep = A*sigtmin*A'+Q;

%% Measurement update
Kt = sigonestep*C'*(C*sigonestep*C'+R)^-1;
meanvarnew.mu = muonestep + Kt*((observedval-d) - C*muonestep);
meanvarnew.sig = sigonestep - Kt*C*sigonestep;
%alpha = 0.5;
%meanvarnew.mu = alpha*observedval+(1-alpha)*meanvar.mu;
%curlocnew = meanvarnew.mu+curlocin;
end
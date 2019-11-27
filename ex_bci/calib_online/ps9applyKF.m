function [meanvarnew] = ps9applyKF(counts,kfparams,meanvar)
A = kfparams.A;
Q = kfparams.Q;
C = kfparams.C;
R = kfparams.R;

if ~exist('meanvar','var')
mu_1_1 = kfparams.pival;
sig_1_1 = kfparams.V;
else
mu_1_1 = meanvar.mu;
sig_1_1 = meanvar.sig;    
end
countsoff = counts - repmat(kfparams.xmean,1,size(counts,2));
%traj = zeros(size(A,1),size(counts,2));
%sigs = cell(1,size(countsoff,2));
%for n = 1:size(countsoff,2)
xt = countsoff; % fill in

% mu_1_1 = mu_0_0; % initialize
% sig_1_1 = sig_0_0; % initialize
mu_1_0 = A*mu_1_1;
sig_1_0 = A*sig_1_1*A'+Q;

kt = sig_1_0*C'*(C*sig_1_0*C'+R)^-1;
meanvarnew.mu = mu_1_0 + kt*(xt - C*mu_1_0);
meanvarnew.sig = sig_1_0 - kt*C*sig_1_0;
% mu_0_0 = (A-kt*C*A)*mu_1_1 + kt*xt
% traj(:,n) = mu_0_0;
% sigs{n} = sig_0_0;
% end
% if exist('z','var')
% mse = mean(sum(((traj(1:2,:) - z(1:2,:))*1000).^2,1))^0.5;
% end
end
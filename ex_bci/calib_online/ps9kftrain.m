function kfparams = ps9kftrain(train_trial,modelparams,sig,sig2,ind1,ind2)
if ~exist('ind1','var'); ind1 = 1; end
if ~exist('ind2','var'); ind2 = 0; end
z1 = [];
z2 = [];
zall = [];
xall = [];
firstzs =[];
for m = 1:length(train_trial)
    if ind2 == 0
        ind2loc = size(train_trial(m).z,2);
    else
        ind2loc = ind2;
    end
    z1 = [z1,train_trial(m).z(:,1:end-1)];
    z2 = [z2,train_trial(m).z(:,2:end)];
    zall = [zall train_trial(m).z(:,ind1:ind2loc)];
    xall = [xall train_trial(m).counts(modelparams.goodneuron==1,ind1:ind2loc)];
    firstzs = [firstzs,train_trial(m).z(:,1)];
end
xmean = mean(xall,2);
xoff = xall-repmat(xmean,1,size(xall,2));
A = (z2*z1')*(z1*z1')^-1;
preq = z2 - A*z1;
if exist('sig2','var')&&~isempty(sig2)
    Q = diag([sig*ones(1,2) sig2*ones(1,2)]);
elseif exist('sig','var')
    Q = diag([sig*ones(1,size(A,1))]);
else
    Q = (1/size(z2,2))*(preq*preq');
end
C = (xoff*zall')*(zall*zall')^-1;
rpre = xoff - C*zall;
R = (1/size(xoff,2))*(rpre*rpre');
pival = mean(firstzs,2);

V = cov(firstzs');
vars = diag(V);
for n=1:length(vars)
    if vars(n)==0
        V(n,n) = 30;
    end
end

kfparams.A = A;
kfparams.Q = Q;
kfparams.C = C;
kfparams.R = R;
kfparams.pival = pival;
kfparams.V = V;
kfparams.xmean = xmean;
end
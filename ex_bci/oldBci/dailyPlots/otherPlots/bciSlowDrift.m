function f=bciSlowDrift(bcistruct,sortflag,bcionlyflag)
if ~exist('bcionlyflag','var');bcionlyflag = 0;end
if ~exist('sortflag','var');sortflag = 0;end
%% for bcistruct
nevcalibtime =cumsum([double(bcistruct.calibrationnev.firstspike);double(bcistruct.calibrationnev.timediff)])/30000;
nevcalib = zeros(length(nevcalibtime),1);
nevcalib(:,1) = double(bcistruct.calibrationnev.nevchan);
nevcalib(:,2) = double(bcistruct.calibrationnev.nevsort);
nevcalib(:,3) = double(nevcalibtime);
nev_infocal = bcistruct.calibrationnev.nev_info;

for n = 1:length(bcistruct.bcinev)
    if n == 1
        nevbcitime = cumsum([(double(bcistruct.bcinev(n).firstspike)+bcistruct.bcinev(n).nev_info.nevclockstart);double(bcistruct.bcinev(n).timediff)])/30000;
        nevbci = zeros(length(nevbcitime),1);
        nevbci(:,1) = double(bcistruct.bcinev(n).nevchan);
        nevbci(:,2) = double(bcistruct.bcinev(n).nevsort);
        nevbci(:,3) = double(nevbcitime);
        nev_infobci(n) = bcistruct.bcinev(n).nev_info;
    else
        nevbcitime = cumsum([(double(bcistruct.bcinev(n).firstspike)+bcistruct.bcinev(n).nev_info.nevclockstart);double(bcistruct.bcinev(n).timediff)])/30000;
        nevbcitemp = zeros(length(nevbcitime),1);
        nevbcitemp(:,1) = double(bcistruct.bcinev(n).nevchan);
        nevbcitemp(:,2) = double(bcistruct.bcinev(n).nevsort);
        nevbcitemp(:,3) = double(nevbcitime);
        nevbci = [nevbci;nevbcitemp];
        nev_infobci(n) = bcistruct.bcinev(n).nev_info;
    end
end

if bcionlyflag == 0
    bignev = [nevcalib;nevbci];
else
    bignev = nevbci;
end
bignevspikes = bignev(bignev(:,1)~=0,:);
if sortflag == 1; bignevspikes = bignevspikes(bignevspikes(:,2)~=0&bignevspikes(:,2)~=255,:);end

stepsize = 60 * 1; % time in seconds
window = 60 * 20; % time in seconds
maxtime = max(bignev(:,3))+stepsize;
mintime = min(bignev(:,3))+stepsize;
[smoothedbins,~,endsteps] = smoothSparse([],bignevspikes,1,stepsize, window,mintime,maxtime);

[coeff, score,latent] = pca(smoothedbins(1:96,:)');

%% plot results
f=figure;
subplot(2,2,1)
plot((endsteps+nev_infocal.nevclockstart/30000)/60,score(:,1))
hold on
endcalibtime = max(nevcalib(:,3));
line([endcalibtime endcalibtime]/60,ylim)
for n = 1:length(nev_infobci)
line([nev_infobci(n).nevclockstart/30000 nev_infobci(n).nevclockstart/30000]/60,ylim)
end
set(gca,'TickDir','out')
box off
xlabel('Time (min)')
ylabel('Projection onto 1st PC')
%figure
subplot(2,2,2)
test = mean(coeff(:,2:end)*score(:,2:end)'+repmat(mean(smoothedbins(1:96,:),2),1,size(smoothedbins(1:96,:),2)),1);
%figure
clear h
h(1)=plot((endsteps+nev_infocal.nevclockstart/30000)/60,mean(smoothedbins(1:96,:),1));
hold on
h(2)=plot((endsteps+nev_infocal.nevclockstart/30000)/60,test);
%ylim([40 60])
line([endcalibtime endcalibtime]/60,ylim)
for n = 1:length(nev_infobci)
line([nev_infobci(n).nevclockstart/30000 nev_infobci(n).nevclockstart/30000]/60,ylim)
end
set(gca,'TickDir','out')
box off
xlabel('Time (min)')
ylabel('Population Firing Rate (Hz)')
legend(h,'With 1st PC','Without 1st PC','Location','Best')
legend boxoff

subplot(2,2,3)
plot(latent/sum(latent))
set(gca,'TickDir','out')
box off
xlabel('PC Index')
ylabel('% Variance Explained')

subplot(2,2,4)
stem(coeff(:,1))
set(gca,'TickDir','out')
box off
xlabel('Coefficient index (PC1)')
ylabel('Coefficient value (PC1)')
end
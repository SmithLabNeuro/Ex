function h =plotSmoothBCIAcquisitionTime(bcidata,maincolor)
if ~exist('maincolor','var');maincolor = 'k';end
bciStarted = [bcidata.bciStarted];
bciStarted_data = bcidata(bciStarted==1);
nolog = zeros(length(bciStarted_data),1);
for n = 1:length(nolog); nolog(n) = isempty(bciStarted_data(n).trialstarttime);end
bciStarted_data = bciStarted_data(nolog==0);

minstepsize = 60 * 1; % time in seconds
window = 60 * 10; % time in seconds
BCI_CORRECT = 161;
bciresults = nan(length(bciStarted_data),1);
trialtime = zeros(length(bciStarted_data),1);
angleval = zeros(length(bciStarted_data),1);
for n = 1:length(bciStarted_data)
    trialcodes = bciStarted_data(n).trialcodes(:,2);
    if isfield(bciStarted_data,'endtrialbin')&&isfield(bciStarted_data,'bciresults')&&bciStarted_data(n).bciresults==1
        bciresults(n) = bciStarted_data(n).endtrialbin;
    elseif ~isempty(find(trialcodes==BCI_CORRECT,1))
        bcistart = bciStarted_data(n).nevinfo.nevclockstart+bciStarted_data(n).trialcodes(find(bciStarted_data(n).trialcodes(:,2)==134,1),3)*30000;
        goodbcibins = bciStarted_data(n).binendtimes>bcistart;
        bciresults(n) = (bciStarted_data(n).binendtimes(end) - bciStarted_data(n).binendtimes((find(goodbcibins,1)-1)))/30000*1000;%add one bcause bcistart is after 1st bin used.
    end
    trialtime(n) = bciStarted_data(n).trialstarttime/30000;
    if isfield(bciStarted_data(n).params.trial,'angle')
        angleval(n) = bciStarted_data(n).params.trial.angle;
    else
        angleval(n) = 0;
    end
end
bcitrials = [bciStarted_data.bcitrial]';
goodresultsinds = ~isnan(bciresults)&bcitrials==1;
goodbciresults = bciresults(goodresultsinds);
goodtrialtimes = trialtime(goodresultsinds);
goodanglevals = angleval(goodresultsinds);
maxtime = max(trialtime);
%mintime = (bcidata(1).nevinfo.nevclockstart/30000-window);
mintime = (bcidata(1).nevinfo.nevclockstart/30000);
[smoothedbins,~,endsteps] = smoothSparse(goodtrialtimes,goodbciresults,0,minstepsize, window,mintime,maxtime);


yaxismax = max(smoothedbins);
h(1) =plot(endsteps/60,smoothedbins,maincolor,'LineWidth',4);
xlabel('time (minutes)')
ylabel(['Acquisition time (ms), window: ',num2str(window/60),' min, step: ',num2str(minstepsize/60),' min'])
hold on
colors = ['rgbm'];
unangles = unique(goodanglevals);
if length(unangles)>1
legendtxt =[{'All Angles'}];
for n = 1:length(unangles)
    thisgoodbciresults = goodbciresults(goodanglevals == unangles(n));
    thisgoodtrialtimes = goodtrialtimes(goodanglevals == unangles(n));
    [smoothedbins,~,endsteps] = smoothSparse(thisgoodtrialtimes,thisgoodbciresults,0,minstepsize, window,mintime,maxtime);
    h(n+1) =plot(endsteps/60,smoothedbins,colors(n),'LineWidth',1);
    thismax=max(smoothedbins);if yaxismax<thismax;yaxismax=thismax;end
    legendtxt{n+1} = ['Angle=',num2str(unangles(n))];
end
legend(h,legendtxt);
%ylim([0 1.05])
legend boxoff
end
set(gca,'TickDir','out')
box off

ylim([bcidata(1).params.block.updatesOnTarget*50 yaxismax+100])
end
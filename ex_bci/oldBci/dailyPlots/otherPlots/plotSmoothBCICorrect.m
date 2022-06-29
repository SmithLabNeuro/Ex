function [h,smoothedbinstruct] = plotSmoothBCICorrect(bcidata,maincolor)
if ~exist('maincolor','var');maincolor = 'k';end
bciStarted = [bcidata.bciStarted];
bciStarted_data = bcidata(bciStarted==1);
nolog = zeros(length(bciStarted_data),1);
for n = 1:length(nolog); nolog(n) = isempty(bciStarted_data(n).trialstarttime);end
bciStarted_data = bciStarted_data(nolog==0);
minstepsize = 60 * 1; % time in seconds
window = 60 * 10; % time in seconds
% aquisition time
% hitrate by time over the session

BCI_CORRECT = 161;
BCI_MISSED = 162;
bciresults = nan(length(bciStarted_data),1);
trialtime = zeros(length(bciStarted_data),1);
angleval = zeros(length(bciStarted_data),1);
for n = 1:length(bciStarted_data)
    if ~isfield(bciStarted_data,'bciresults')
        trialcodes = bciStarted_data(n).trialcodes(:,2);
        if ~isempty(find(trialcodes==BCI_CORRECT,1))
            bciresults(n) = 1;
        elseif ~isempty(find(trialcodes==BCI_MISSED,1))
            bciresults(n) = 0;
        end
    else
        bciresults(n) = bciStarted_data(n).bciresults;
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

overallsmoothed = smoothedbins;
h(1) =plot(endsteps/60,smoothedbins,maincolor,'LineWidth',4);
xlabel('time (minutes)')
ylabel(['% Correct, window: ',num2str(window/60),' min, step: ',num2str(minstepsize/60),' min'])
hold on
colors = ['rgbm'];
unangles = unique(goodanglevals);
if length(unangles)>1
legendstrings = [{'All Angles'}];
for n = 1:length(unangles)
    thisgoodbciresults = goodbciresults(goodanglevals == unangles(n));
    thisgoodtrialtimes = goodtrialtimes(goodanglevals == unangles(n));
    [smoothedbins,~,endsteps] = smoothSparse(thisgoodtrialtimes,thisgoodbciresults,0,minstepsize, window,mintime,maxtime);
    h(n+1) =plot(endsteps/60,smoothedbins,colors(n),'LineWidth',1);
    legendstrings{n+1} = ['Angle=',num2str(unangles(n))];
    smoothedbinstruct(n).angle = unangles(n);
    smoothedbinstruct(n).smoothedbins = smoothedbins;
    smoothedbinstruct(n).unsmoothedbins = thisgoodbciresults;
    smoothedbinstruct(n).thisgoodtrialtimes = thisgoodtrialtimes;    
end

smoothedbinstruct(end+1).smoothedbins = overallsmoothed;
legend(h,legendstrings,'Location','Best')
legend boxoff
end
set(gca,'TickDir','out')
box off

ylim([0 1.05])
end

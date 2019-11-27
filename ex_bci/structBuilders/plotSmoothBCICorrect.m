function h = plotSmoothBCICorrect(bcidata)
bciStarted = [bcidata.bciStarted];
bciStarted_data = bcidata(bciStarted==1);

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
    angleval(n) = bciStarted_data(n).params.trial.angle;
end
bcitrials = [bciStarted_data.bcitrial]';
goodresultsinds = ~isnan(bciresults)&bcitrials==1;
goodbciresults = bciresults(goodresultsinds);
goodtrialtimes = trialtime(goodresultsinds);
goodanglevals = angleval(goodresultsinds);
maxtime = max(trialtime);
mintime = (bcidata(1).nevinfo.nevclockstart/30000-window);
[smoothedbins,~,endsteps] = smoothSparse(goodtrialtimes,goodbciresults,0,minstepsize, window,mintime,maxtime);


h(1) =plot(endsteps/60,smoothedbins,'k','LineWidth',4);
xlabel('time (minutes)')
ylabel(['% Correct, window: ',num2str(window/60),' min, step: ',num2str(minstepsize/60),' min'])
hold on
colors = ['rgbm'];
unangles = unique(goodanglevals);
for n = 1:4
    thisgoodbciresults = goodbciresults(goodanglevals == unangles(n));
    thisgoodtrialtimes = goodtrialtimes(goodanglevals == unangles(n));
    [smoothedbins,~,endsteps] = smoothSparse(thisgoodtrialtimes,thisgoodbciresults,0,minstepsize, window,mintime,maxtime);
    h(n+1) =plot(endsteps/60,smoothedbins,colors(n),'LineWidth',1);
end
set(gca,'TickDir','out')
box off
legend(h,'All Angles',['Angle=',num2str(unangles(1))],['Angle=',num2str(unangles(2))],['Angle=',num2str(unangles(3))],['Angle=',num2str(unangles(4))],'Location','Best')
legend boxoff
ylim([0 1.05])
end

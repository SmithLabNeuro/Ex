function otherDailyPlots(bcistruct,filename,newsavepath)
if ~exist('newsavepath','var')
    newsavepath = [];
end
    tmp = strsplit(filename,'_');
    subjectDate = tmp{1};
figure
subplot(2,1,1)
plotSmoothBCICorrect(bcistruct.bcidata);
subplot(2,1,2)
plotSmoothBCIAcquisitionTime(bcistruct.bcidata);
set(gcf,'Position',[100 100 800 800])
[~,name] = fileparts(filename);
saveas(gcf,[newsavepath,'/',name,'performanceOverTime'],'png')
% f=bciSlowDrift(bcistruct,1,1);
% saveas(f,[newsavepath,'/','nevPCASlowDrift'],'png')
end
function otherDailyPlots(bcistruct,filename)
    tmp = strsplit(filename,'_');
    subjectDate = tmp{1};
figure
subplot(2,1,1)
plotSmoothBCICorrect(bcistruct.bcidata);
subplot(2,1,2)
plotSmoothBCIAcquisitionTime(bcistruct.bcidata);
set(gcf,'Position',[100 100 800 800])
saveas(gcf,[subjectDate,'/','performanceOverTime'],'png')
end
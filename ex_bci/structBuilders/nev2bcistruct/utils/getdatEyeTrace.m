function [eyedeg,eyepix] = getdatEyeTrace(dat,trialnum,plotflag)
 eyedeg=eye2deg(dat(trialnum).eyedata.trial,dat(trialnum).params);
 eyepix = deg2pix(eyedeg,dat(trialnum).params.block.screenDistance,dat(trialnum).params.block.pixPerCM);
 
 if plotflag == 1
    figure;plot(eyepix(1,:),eyepix(2,:))
    xlim([-dat(trialnum).params.block.displayWidth/2 dat(trialnum).params.block.displayWidth/2])
    ylim([-dat(trialnum).params.block.displayHeight/2 dat(trialnum).params.block.displayHeight/2])
 end
end
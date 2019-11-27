function [dat] = dat2Traj(dat,modelparams)
if ~exist('modelparams','var')
    recalibmodelflag = 1;
else
    recalibmodelflag = 0;
end
startval = 1;
%%%% outside code for recalibrating bci %%%%

    numcalibtrials = 800;
    recalibinterval = 2;
if recalibmodelflag == 1
    startval = numcalibtrials+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targsize = 80;
%% generate trajectories
for n = startval:length(dat)
    if mod(n,100)==0
        fprintf('processing trial %d\n',n)
    end
    spikecounts = dat(n).counts;
    %%%%%%%% Code for recalibrating bci %%%%%%%%%%%
    recalibflag =mod(n-numcalibtrials,recalibinterval);
    if recalibmodelflag == 1 && recalibflag == 1
        countmat = [];
        anglemat = [];
        for m = (n-numcalibtrials):(n-1)
            countmat = [countmat dat(m).counts];
            anglemat = [anglemat dat(m).angle*ones(1,size(dat(m).counts,2))];
        end    
        uniqueangle = unique(anglemat);
        modelparams.angles = uniqueangle;
        modelparams.step = 13;%0.1*max(distancemat);
        meancounts = zeros(size(countmat,1),length(uniqueangle));
        for m = 1:length(uniqueangle)
            meancounts(:,m) = mean(countmat(:,anglemat==uniqueangle(m)),2);
        end
        modelparams.meancounts = meancounts;
        modelparams.goodneuron = mean(countmat,2)>(1*0.1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cursorpos1 = 0;
    cursorpos2 = 0;
    %meanvar.mu = 0; meanvar.sig = 1;
    dat(n).angleerror = 0;
    %meanvar2.mu = 0; meanvar2.sig = 1;
    for m = 1:size(spikecounts,2)
    %[cursorpos1temp, cursorpos2temp] = decodeNBPCursor(spikecounts(modelparams.goodneuron==1,m),modelparams);
    if m == 1
        [meanvar ] = ps9applyKF(spikecounts(modelparams.goodneuron==1,m),modelparams.KFparamsdec);
    else
        [meanvar ] = ps9applyKF(spikecounts(modelparams.goodneuron==1,m),modelparams.KFparamsdec,meanvar);
        modelparams.KFparamsdec.V = meanvar.sig;
    end
    cursorpos1(m) = meanvar.mu(1);
    cursorpos2(m) = meanvar.mu(2);
    theta = deg2rad(dat(n).angle);
    dat(n).targX = round(130*cos(theta));
    dat(n).targY = round(130*sin(theta));
    %% add code for computing error
    muvec = [meanvar.mu(1) meanvar.mu(2)];
    targvec = [dat(n).targX dat(n).targY];
    dat(n).angleerror(m)= rad2deg(acos(dot(muvec,targvec)/(norm(muvec)*norm(targvec))));
    
    end
    %%%%%%%%%%%%%% Set target as fixed (constant distance across trials %%%
     theta = deg2rad(dat(n).angle);
    dat(n).targX = round(130*cos(theta));
    dat(n).targY = round(130*sin(theta));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dat(n).traj1 = [0 cumsum(cursorpos1)];%[0 cursorpos1];%
    dat(n).traj2 = [0 cumsum(cursorpos2)];%[0 cursorpos2];%
    dat(n).mse = mean(sum((([dat(n).traj1; dat(n).traj2] - [[0;0] cumsum(dat(n).z(1:2,:),2)])*1000).^2,1))^0.5;
    dat(n).meananglerr = mean(dat(n).angleerror);
    if sum((((dat(n).traj1-dat(n).targX).^2 + (dat(n).traj2-dat(n).targY).^2).^0.5)<targsize)>0   
        dat(n).bcicorrect = 1;
    else
        dat(n).bcicorrect = 0;
    end
end
%% find condition value
if recalibmodelflag == 1
    datplot = dat(801:end);
else
    datplot = dat;
end
angle = [datplot.angle];
unangle = unique(angle);


%% plot trajectory by condition value
figure
subplotrows = 3;
subplotcol = 3;
subplotinds = [6 3 2 1 4 7 8 9];%[1 2 3 4 6 7 8 9];
a=distinguishable_colors(length(unangle)+1);
plotsize = 150;
for n = 1:length(unangle)
    tempdat = datplot(angle==unangle(n));    
    for m = 1:length(tempdat)

        subplot(subplotrows,subplotcol,subplotinds(n))
        plot(tempdat(m).traj1,tempdat(m).traj2,'Color',a(n,:));
        hold on
        xlim([-1*plotsize plotsize])
        ylim([-1*plotsize plotsize])
        subplot(subplotrows,subplotcol,5)
        plot(tempdat(m).traj1,tempdat(m).traj2,'Color',a(n,:));
        hold on
        xlim([-1*plotsize plotsize])
        ylim([-1*plotsize plotsize])
        
    end
    targlocx = tempdat(1).targX;
    targlocy = tempdat(1).targY;
    subplot(subplotrows,subplotcol,subplotinds(n))
    t = linspace(0,2*pi);plot(targsize*cos(t)+targlocx,targsize*sin(t)+targlocy,'Color',a(end,:),'LineWidth',2)
end



end
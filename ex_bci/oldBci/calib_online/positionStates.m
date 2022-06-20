function [Xmat,Ymat,dat] = positionStates(dat,numsteps,multiplier,reallengthflag)
targx = [dat.newX]*multiplier;
targy = [dat.newY]*multiplier;
Xmat = [];
Ymat = [];
for n = 1:length(targx)
    if reallengthflag == 0
        x = linspace(0,targx(n),numsteps);
        y = linspace(0,targy(n),numsteps);
        Xmat = [Xmat x];
        Ymat = [Ymat y];
        dat(n).z = [ones(1,length(x))*(x(2)-x(1));ones(1,length(y))*(y(2)-y(1))]; 
    else
        x = linspace(0,targx(n),numsteps);
        y = linspace(0,targy(n),numsteps);
        Xmat = [Xmat x];
        Ymat = [Ymat y];
        dat(n).z = [ones(1,size(dat(n).counts,2))*(x(2)-x(1));ones(1,size(dat(n).counts,2))*(y(2)-y(1))]; 
    end
end

end


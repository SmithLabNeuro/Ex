function outinds = binaryTimeSearch(nevtime,tstart,tend)
% tend is optional. If not input, then will only output the index of tstart

converged = 0;
maxind = length(nevtime);
minind = 1;
startind = round(maxind/2);
while converged == 0
    if nevtime(startind) < tstart
        minind = startind;
        startind = round((maxind-minind)/2) + minind;
    else
        maxind = startind;
        startind = maxind-round((maxind-minind)/2);
    end
    if (maxind-minind) < 2
        converged = 1;
    end
end
finalinds =[minind maxind];
finalstart = finalinds(find(nevtime(finalinds)>=tstart,1));
%finalstart = binarySearch(nevtime,tstart,'first','up');
%finalend = binarySearch(nevtime,tend,'first','down');
%outinds = (finalstart:finalend);
if exist('tend','var')
converged = 0;
maxind = length(nevtime);
minind = 1;
startind = round(maxind/2);
while converged == 0
    if nevtime(startind) < tend
        minind = startind;
        startind = round((maxind-minind)/2) + minind;
    else
        maxind = startind;
        startind = maxind-round((maxind-minind)/2);
    end
    if (maxind-minind) < 2
        converged = 1;
    end
end
finalinds =[minind maxind];
finalend = finalinds(find(nevtime(finalinds)<tend,1,'last'));
outinds = (finalstart:finalend);
else
    outinds = finalstart;
end
end
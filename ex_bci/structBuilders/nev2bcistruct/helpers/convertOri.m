function trueori = convertOri(oriPick, isValid, cue)
    trueori = zeros(length(oriPick),1);
    for n = 1:length(oriPick)
        tempOri = oriPick(n); %'oriPick' is the XML parameter that is sent for the trial (refers to target stimulus)
        if cue(n)==2, %if cue is out of RF...
            tempOri = 3-tempOri; %...flip the ori to refer to the RF.
        end;
        if isValid(n), %if the cue is invalid...
            tempOri = 3-tempOri; %...flip the ori to refer to the uncued location.
        end;
        trueori(n) = tempOri; %this is now the orientation at the RF
    end
end
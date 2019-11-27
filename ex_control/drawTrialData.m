function drawTrialData
% function drawTrialData
%
% simply writes the current trialData to the screen. used by runex and can
% be used within an ex file

global trialData wins params statusInfo;

x = 0;
y = 0;
gray=(WhiteIndex(0)+BlackIndex(0))/2;

Screen(wins.info,'FillRect',gray);
Screen('TextSize',wins.info,wins.textSize);
for i = 1:length(trialData)
    if ~isempty(trialData{i}),
        Screen('DrawText',wins.info,trialData{i},x,y);
    end;
    y = y + wins.lineSpacing.*wins.textSize;
end

%Write a small file to the status folder for the monitoring applet
%(thanks, Mike Morais!) -06Sep2013 ACS:
if params.statusUpdates
    try
        aioV = samp;
        aioV = aioV(end,:);
        eyeLine = sprintf('eyeX: %.2f, eyeY: %.2f',aioV);
        fid = fopen(statusInfo.filename,'w');
        if fid>-1
            fprintf(fid,'%s\n',datestr(now),params.SubjectID,eyeLine,trialData{:});
            fclose(fid);
        end;
    catch ME
        warning('Error updating status:: %s',ME.message); %#ok<WNTAG>
    end;
end;

end

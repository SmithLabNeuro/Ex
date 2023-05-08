function drawTrialData
% function drawTrialData
%
% simply writes the current trialData to the screen. used by runex and can
% be used within an ex file

global trialData wins params;

x = 0;
y = 0;
gray=(WhiteIndex(0)+BlackIndex(0))/2;

Screen(wins.info,'FillRect',gray);
Screen('TextSize',wins.info,wins.textSize);
for i = 1:length(trialData)
    if ~isempty(trialData{i})
        if size(wins.trialData.lineColor,1)>= i
            textColor = wins.trialData.lineColor(i,:);
        else
            textColor = [0,0,0]; % default to black for extra lines
        end
        Screen('DrawText',wins.info,trialData{i},x,y,textColor);
    end
    y = y + wins.lineSpacing.*wins.textSize;
end

end

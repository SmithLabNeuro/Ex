function [success, msg] = bciUpdatePhotoDiode(~,~, bciSocket, photoDiodeObj, photoDiodeSquarePosition, photoDiodeSquareWidth)

global behav

if matlabUDP2('check', bciSocket)
    receivedMsg = matlabUDP2('receive', bciSocket);
    if ~isempty(receivedMsg)
        %         disp(receivedMsg)
        meanSpikeCount = typecast(uint8(receivedMsg), 'double');
        maxMeanCount = 10;%x0.5;
        minMeanCount = 0;
        if meanSpikeCount > maxMeanCount
            meanSpikeCount = maxMeanCount;
        end
        
        % define the intensity from 0-255--this will be the RGB for a grayscale
        % color
        photoDiodeIntensity = round(255*(meanSpikeCount - minMeanCount)/(maxMeanCount - minMeanCount));
        if isnan(photoDiodeIntensity)
%             disp(int8(receivedMsg))
%             disp(uint8(receivedMsg))
            disp('nan')
            %             keyboard
            msg = '';
        else
%             disp(int8(receivedMsg))
%             disp(uint8(receivedMsg))
            msg = sprintf('set %i rect 0 %i %i %i %i %i %i %i', photoDiodeObj, photoDiodeSquarePosition(1), photoDiodeSquarePosition(2), photoDiodeSquareWidth, photoDiodeSquareWidth, photoDiodeIntensity, photoDiodeIntensity, photoDiodeIntensity);
        end
        if ~isfield(behav, 'pdIntensityMeanCount')
            behav.pdIntensityMeanCount = nan(100000, 1);
        end
        behav.pdIntensityMeanCount = meanSpikeCount;
    else
        disp('emp')
        msg = '';
    end
    % this function is always successful--it's a nonblocking function on
    % the loop
else
    msg = '';
end
success = 1;
end
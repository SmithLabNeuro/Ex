function stim_movie(optstr,w,objID,arg)
%function stim_movie(optstr,w,objID,arg)
%
% showex helper function for 'movie' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    [fileName, argsRest] = strtok(arg);
    a = sscanf(argsRest,'%i');
    % arguments: (1) frameCount
    %            (2) dwell
    %            (3) start frame
    %            (4) x position
    %            (5) y position
    
    % additional options: (6) width     %added 20dec2018
    %                     (7) height
    
    dwell = a(2);
    startFrame = a(3);
    xCenter = a(4);
    yCenter = a(5) * -1; % flip y coordinate so '-' is down
    
    
    
    vars = load([sv.localDir,filesep,fileName]); %added sv.localDir, now prefix/suffix combos will be relative to that location. This directory can be set by rig/subject XMLs. -ACS 13Sep2013
    
    %adding 'dstSize' so that the movie can be resized by PTB
    %-acs20dec2018:
    if numel(a) > 5
        if numel(a) == 6
            a(7) = a(6);
        end;
        dstSize = a(6:7)';
    else
        dstSize = size(vars.mov{1});
    end;
    movieSize = size(vars.mov{1});
    % MAS 06/27/13 - removed the "-1" from
    % dstRect/srcRect. Also made the code use the first
    % two values in movieSize for when a 3D (color)
    % movie is in use
    dstRect=[0 0 dstSize(2:-1:1)];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];
    srcRect=[0 0 movieSize(2:-1:1)];

    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1).*dwell, ...    %changed so that the total frame count is the number of movie frames times the dwell -acs20dec2018
        'startFrame',startFrame, 'dwell',dwell, ...
        'srcRect',srcRect, 'dstRect',dstRect);
    
    objects{objID}.mov = vars.mov;
elseif strcmp(optstr,'display')
    frameNum = floor(objects{objID}.frame/objects{objID}.dwell) + objects{objID}.startFrame;
    
    if frameNum > length(objects{objID}.mov)
        frameNum = mod(frameNum-1,length(objects{objID}.mov))+1;
    end

    tex=Screen('MakeTexture', w, objects{objID}.mov{frameNum});
    Screen('DrawTexture', w, tex,objects{objID}.srcRect,objects{objID}.dstRect);
    %                    Screen('DrawTexture', arg, tex);
    
    Screen('Close',tex);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

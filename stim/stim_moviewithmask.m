function stim_moviewithmask(optstr,w,objID,arg)
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
    [fileName argsRest] = strtok(arg);
    a = sscanf(argsRest,'%i %i %i %i %i %f %i');
    
    % arguments: (1) frameCount
    %            (2) dwell
    %            (3) start frame
    %            (4) x position
    %            (5) y position
    
    dwell = a(2);
    startFrame = a(3);
    xCenter = a(4);
    yCenter = a(5) * -1; % flip y coordinate so '-' is down
    contrast = a(6);
    maskFrame = a(7);
    vars = load([sv.localDir,filesep,fileName]); %added sv.localDir, now prefix/suffix combos will be relative to that location. This directory can be set by rig/subject XMLs. -ACS 13Sep2013
    
    movieSize = size(vars.mov{1});
    % MAS 06/27/13 - removed the "-1" from
    % dstRect/srcRect. Also made the code use the first
    % two values in movieSize for when a 3D (color)
    % movie is in use
    dstRect=[0 0 movieSize(2:-1:1)]
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter]
    srcRect=[0 0 movieSize(2:-1:1)]
    randframe = randi(350);
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'startFrame',startFrame, 'dwell',dwell, ...
        'srcRect',srcRect, 'dstRect',dstRect,'contrast',contrast,'maskFrame',maskFrame,'randFrame',randframe);
    
    objects{objID}.mov = vars.mov;
    
elseif strcmp(optstr,'display')
    frameNum = floor(objects{objID}.frame/objects{objID}.dwell) + objects{objID}.startFrame;
    
    if frameNum > length(objects{objID}.mov)
        frameNum = mod(frameNum-1,length(objects{objID}.mov))+1;
    end
    objects{objID}.maskFrame
    frameNum
    randframe = objects{objID}.randFrame;
    %tex=Screen('MakeTexture', w, objects{objID}.mov{frameNum});   
    tex=Screen('MakeTexture', w, objects{objID}.mov{randframe}); 
    Screen('DrawTexture', w, tex,objects{objID}.srcRect,objects{objID}.dstRect);
    Screen('Close',tex);
    if objects{objID}.maskFrame<=frameNum && objects{objID}.contrast > 0
       % [s1,s2,s3] = size(objects{objID}.mov{frameNum});
        [s1,s2,s3] = size(objects{objID}.mov{randframe});
        mask = ones(s1, s2, 2) .* sv.bgColor(1); % make this the background color (maybe make this an input)
        mask(:, :, 2)= round(255*ones(s1,s2,1).*objects{objID}.contrast); % contrast can be an input?
        maskTexture = Screen('MakeTexture', w, mask);
        Screen('DrawTexture', w, maskTexture,objects{objID}.srcRect,objects{objID}.dstRect);
        Screen('Close',maskTexture);
    end
    
    
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end
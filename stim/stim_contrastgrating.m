function stim_contrastgrating(optstr,w,objID,arg)
%
% showex helper function for 'grating' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %f %f %i %i %i %i %i %f %i %f %f %i');
    
    % arguments: (1) frameCount
    %            (2) angle
    %            (3) spatial frequency
    %            (4) x position
    %            (5) y position
    %            (6) aperture size
    %            (7) num flashes before change
    %            (8) frames per contrast cycle
    %            (9) amplitdue of sine contrast
    %            (10) bool, use phase swap?
    %            (11) base hue H
    %            (12) new hue H
    %            (13) hues alpha
    a
    angle = mod(180-a(2),360);
    f = a(3);
    xCenter = a(4);
    yCenter = -a(5); % flip y coordinate so '-' is down
    rad = a(6); % Size of the grating image. Needs to be a power of two.
    
    % Calculate parameters of the grating:
    ppc=ceil(1/f);  % pixels/cycle
    fr=f*2*pi; % radians
    visibleSize=2*rad+1;
    
    % Create one single static grating image:
    x=meshgrid(-rad:rad + ppc, -rad:rad);
    grating = sv.gray + (sv.gray*cos(fr*x));
    
    % Store grating in texture: Set the 'enforcepot' flag to 1 to signal
    % Psychtoolbox that we want a special scrollable power-of-two texture:
    gratingTex=Screen('MakeTexture', w, grating);
    
    % Create a single gaussian transparency mask and store it to a texture:
    % this makes the visible grating a circle and not a square.
    mask=ones(2*rad+1, 2*rad+1, 2) * mean(sv.bgColor);
    [x,y]=meshgrid(-1*rad:1*rad,-1*rad:1*rad);
    
    mask(:, :, 2)=sv.white * (sqrt(x.^2+y.^2) > rad);
    
    maskTex=Screen('MakeTexture', w, mask);
    
    % contrast modulation sine wave parameters
    amplitude = a(9);
    contrastStepRad = pi / a(8);
    startPhase = 0;
    numFlashes = a(7);
    
    % color mask parameters
    saturation = 1;
    value = 0.67;
    initialHsv = [a(11), saturation, value];
    changedHsv = [a(12), saturation, value];
    initialRgb = round(hsv2rgb(initialHsv)*255);
    changedRgb = round(hsv2rgb(changedHsv)*255);
    initialCol = [initialRgb, a(13)]; % RGBA
    changedCol = [changedRgb, a(13)];

    dstRect=[0 0 visibleSize visibleSize];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'angle',angle,'size',visibleSize,'x',xCenter,'y',yCenter,'phaseOffset',[0, ceil(ppc/2)],'phaseIdx',1,'usePhaseSwap',logical(a(10)),...
        'sineAmp',amplitude,'sineStep',contrastStepRad,'sinePhase',startPhase,'numFlashes',numFlashes,'numSteps',a(8),...
        'grating',gratingTex, 'mask',maskTex, ...
        'ppc',ppc,'dstRect',dstRect,'initialCol',initialCol,'changedCol',changedCol,'rad',rad);
elseif strcmp(optstr,'display')
    % check if contrast is 0 and phase needs to be swapped
    if objects{objID}.usePhaseSwap
        if mod(objects{objID}.frame, objects{objID}.numSteps) == 0
            tmp = objects{objID}.frame / objects{objID}.numSteps;
            if mod(tmp,2) == 1
                objects{objID}.phaseIdx = mod(objects{objID}.phaseIdx,2) + 1;
            end
        end
    end
    % check if color needs to be changed
    if objects{objID}.frame >= (objects{objID}.numSteps) * (2 * objects{objID}.numFlashes - 1)
        currCol = objects{objID}.changedCol;
    else
        currCol = objects{objID}.initialCol;
    end
    currContrast = objects{objID}.sineAmp * cos(objects{objID}.sineStep * objects{objID}.frame + objects{objID}.sinePhase) + (1 - objects{objID}.sineAmp);

    xOffset = objects{objID}.phaseOffset(objects{objID}.phaseIdx);
    srcRect = [xOffset 0 xOffset + objects{objID}.size objects{objID}.size];
    targetPos = [sv.midScreen + [objects{objID}.x objects{objID}.y] - objects{objID}.rad,sv.midScreen + [objects{objID}.x objects{objID}.y] + objects{objID}.rad];

    Screen('DrawTexture',w,objects{objID}.grating,srcRect,objects{objID}.dstRect,objects{objID}.angle, [], currContrast);
    Screen('DrawTexture',w,objects{objID}.mask,[0 0 objects{objID}.size objects{objID}.size],objects{objID}.dstRect,objects{objID}.angle);
    Screen(w,'FillOval',currCol,targetPos);
elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.grating);
    Screen('Close',objects{objID}.mask);
else
    error('Invalid option string passed into stim_*.m function');
end

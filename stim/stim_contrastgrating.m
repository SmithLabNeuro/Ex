function stim_grating(optstr,w,objID,arg)
%function stim_grating(optstr,w,objID,arg)
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
    a = sscanf(arg,'%i %f %f %f %f %i %i %i %f %i %f %f');
    
    % arguments: (1) frameCount
    %            (2) angle
    %            (3) initial phase
    %            (4) frequency
    %            (5) cycles per second
    %            (6) x position
    %            (7) y position
    %            (8) aperture size
    %            (9) contrast (0.0-1.0)
    %            (10) num flashes before change
    %            (11) frames per contrast cycle
    %            (12) amplitdue of sine contrast
    
    angle = mod(180-a(2),360);
    mult_angle = [0, angle];
    f = a(4);
    cps = a(5);
    xCenter = a(6);
    yCenter = -a(7); % flip y coordinate so '-' is down
    rad = a(8); % Size of the grating image. Needs to be a power of two.
    contrast = a(9);
    
    % Calculate parameters of the grating:
    ppc=ceil(1/f);  % pixels/cycle
    fr=f*2*pi; % radians
    visibleSize=2*rad+1;
    
    phase = a(3)/360*ppc;
    
    % Create one single static grating image:
    x=meshgrid(-rad:rad + ppc, -rad:rad);
    grating = sv.gray + (sv.gray*cos(fr*x)) * contrast;
    
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
    amplitude = a(12);
    contrastStepRad = pi / a(11);
%     contrastFreqRad = contrastFreq * 2 * pi;
    startPhase = 0;
    numFlashes = a(10);
    
    dstRect=[0 0 visibleSize visibleSize];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'angles',mult_angle, 'phase',phase, ...
        'size',visibleSize, 'x',xCenter,'y',yCenter,...
        'sineAmp', amplitude, 'sineStep', contrastStepRad, 'sinePhase', startPhase, 'numFlashes', numFlashes, 'numSteps', a(11),...
        'grating',gratingTex, 'mask',maskTex, ...
        'ppc',ppc, 'dstRect',dstRect);
elseif strcmp(optstr,'display')
    if objects{objID}.frame >= (objects{objID}.numSteps) * (2 * objects{objID}.numFlashes - 1)
        currAngle = objects{objID}.angles(2);
    else
        currAngle = objects{objID}.angles(1);
    end
    currContrast = objects{objID}.sineAmp * cos(objects{objID}.sineStep * objects{objID}.frame + objects{objID}.sinePhase) + (1 - objects{objID}.sineAmp);
    srcRect = [0 0 objects{objID}.size objects{objID}.size];
    
    Screen('DrawTexture',w,objects{objID}.grating,srcRect,objects{objID}.dstRect,currAngle, [], currContrast);
    Screen('DrawTexture',w,objects{objID}.mask,[0 0 objects{objID}.size objects{objID}.size],objects{objID}.dstRect,currAngle);
elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.grating);
    Screen('Close',objects{objID}.mask);
else
    error('Invalid option string passed into stim_*.m function');
end

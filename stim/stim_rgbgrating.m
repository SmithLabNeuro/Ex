function stim_rgbgrating(optstr,w,objID,arg)
%function stim_rgbgrating(optstr,w,objID,arg)
%
% showex helper function for 'rgbgrating' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %f %f %f %f %i %i %i %f %i %i %i %i %i %i %i %f %i %i %i');
    
    % arguments: (1) frameCount
    %            (2) angle
    %            (3) initial phase
    %            (4) frequency
    %            (5) cycles per second (baseline)
    %            (6) x position
    %            (7) y position
    %            (8) aperture size
    %            (9) contrast (0.0-1.0)
    %           (10) start colormap R
    %           (11) end colormap R
    %           (12) start colormap G
    %           (13) end colormap G
    %           (14) start colormap B
    %           (15) end colormap B
    %           (16) colormap resolution
    %           (17) peak cps
    %           (18) start ramp frame
    %           (19) ramp peak frame
    %           (20) end ramp frame
    
    frameCount = a(1);
    angle = mod(180-a(2),360);
    f = a(4);
    baseCps = a(5); %-ACS20Feb2012
    peakCps = a(17); %-ACS20Feb2012
    xCenter = a(6);
    yCenter = -a(7); % flip y coordinate so '-' is down
    rad= a(8); % Size of the grating image. Needs to be a power of two.
    contrast = a(9);
    startRamp = a(18);
    peakRamp = a(19);
    endRamp = a(20);
    
    % Calculate parameters of the grating:
    ppc=ceil(1/f);  % pixels/cycle
    fr=f*2*pi;
    visibleSize=2*rad+1;
    
    phase = a(3)/360*ppc;
    
    % Create one single static grating image:
    x=meshgrid(-rad:rad + ppc, -rad:rad);
    grating = sv.gray + (sv.inc*cos(fr*x))*contrast;
    
    % Impose colormap (ACS20Feb2012):
    cmap = [linspace(a(10),a(11),a(16))',linspace(a(12),a(13),a(16))',linspace(a(14),a(15),a(16))'];
    grating = ind2rgb(gray2ind(grating./max(grating(:)),a(16)),cmap);
    
    % Store grating in texture: Set the 'enforcepot' flag to 1 to signal
    % Psychtoolbox that we want a special scrollable power-of-two texture:
    gratingTex=Screen('MakeTexture', w, grating);
    
    % Create a single gaussian transparency mask and store it to a texture:
    mask=ones(2*rad+1, 2*rad+1, 2) * mean(sv.bgColor);
    [x,y]=meshgrid(-1*rad:1*rad,-1*rad:1*rad);
    
    mask(:, :, 2)=sv.white * (sqrt(x.^2+y.^2) > rad);
    
    maskTex=Screen('MakeTexture', w, mask);
    
    %Make shift function (ACS20Feb2011):
    shiftFun = baseCps*ones(1,frameCount);
    shiftFun(startRamp:peakRamp-1) = linspace(baseCps,peakCps,peakRamp-startRamp);
    shiftFun(peakRamp:endRamp-1) = linspace(peakCps,baseCps,endRamp-peakRamp);
    shift = cumsum(shiftFun.*ppc.*sv.ifi);
    
    dstRect=[0 0 visibleSize visibleSize];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];

    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'angle',angle, 'phase',phase, 'shift', shift, ...
        'size',visibleSize, 'x',xCenter,'y',yCenter, ...
        'grating',gratingTex, 'mask',maskTex, ...
        'ppc',ppc, 'dstRect',dstRect);
elseif strcmp(optstr,'display')
    if isscalar(objects{objID}.shift) %added dynamic shift support (ACS20Feb2011)
        xOffset = mod(objects{objID}.frame*objects{objID}.shift+objects{objID}.phase,objects{objID}.ppc);
    else
        xOffset = mod(objects{objID}.shift(objects{objID}.frame+1)+objects{objID}.phase,objects{objID}.ppc); %don't multiply by frame b/c of cumsum in shiftFun -ACS20Feb2012
    end;
    srcRect = [xOffset 0 xOffset + objects{objID}.size objects{objID}.size];
    
    Screen('DrawTexture',w,objects{objID}.grating,srcRect,objects{objID}.dstRect,objects{objID}.angle);
    Screen('DrawTexture',w,objects{objID}.mask,[0 0 objects{objID}.size objects{objID}.size],objects{objID}.dstRect,objects{objID}.angle);
elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.grating);
    Screen('Close',objects{objID}.mask);
else
    error('Invalid option string passed into stim_*.m function');
end

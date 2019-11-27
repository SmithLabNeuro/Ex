function stim_gabor(optstr,w,objID,arg)
%function stim_gabor(optstr,w,objID,arg)
% NOTE: This method for displaying gabors doesn't work as of 01oct2015
% because our Linux version is too out of date and is missing needed
% libraries for OpenGL. -ACS 01oct2015
%
% showex helper function for 'gabor' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %f %f %f %f %i %i %i %f %f %i %i %i');
    
    if numel(a)==10, %indicates that no bg color was passed - for backwards compatibility -acs22Jun2016
        a = [a(:)' 128 128 128];  %assume bg color is 50% gray -acs22jun2016  
    end;
    
    % arguments: (1) frameCount
    %            (2) angle
    %            (3) initial phase
    %            (4) frequency
    %            (5) cycles per second
    %            (6) x position
    %            (7) y position
    %            (8) radius
    %            (9) contrast (0.0-1.0)
    %           (10) sigma of gaussian
    %           (11) bgColor 1
    %           (12) bgColor 2
    %           (13) bgColor 3
    
    a(2) = -a(2); %Gabors rotate in the opposite directionfrom gratings with increasing angles, so negate the angle here for consistency. -acs17jan2017
    
    angle = mod(180-a(2),360);
    f = a(4);
    cps = a(5);
    xCenter = a(6);
    yCenter = -a(7); % flip y coordinate so '-' is down
    rad= a(8); % Size of the grating image. Needs to be a power of two.
    contrast = a(9); %rescale to percentages? -acs01oct2015
    sc = a(10);
    
    % Calculate parameters of the grating:
    visibleSize=2*rad+1;
    
    phase = a(3);
    
    shift = 360 * cps * sv.ifi;
    
    dstRect=[0 0 visibleSize visibleSize];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];
    
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'angle',angle, 'phase',phase, 'shift',shift, ...
        'size',visibleSize, 'x',xCenter,'y',yCenter, ...
        'freq',f,'sc',sc,'contrast',contrast, 'dstRect',dstRect);
    
    objects{objID}.gabortex = CreateProceduralGabor(w, objects{objID}.size, objects{objID}.size, 0, [a(11)/255 a(12)/255 a(13)/255 1],1,0.5); %nb: will need to fix the background color to be detected... going with 0.5's for now... -acs 01oct2015

elseif strcmp(optstr,'display'),    
    thisPhase = objects{objID}.frame*objects{objID}.shift+objects{objID}.phase;  
    Screen('DrawTexture',w,objects{objID}.gabortex,[],objects{objID}.dstRect,objects{objID}.angle,[],[],[255,255,255,0],[], kPsychDontDoRotation, [thisPhase, objects{objID}.freq, objects{objID}.sc, objects{objID}.contrast, 1, 0, 0, 0]); %The "modulateColor" parameter ranges from 0 (no modulation) to 255 (full modulation). Note that the alpha channel is zero here, and a constant is added to the base color when the gabor is defined (above). -acs01oct2015

elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.gabortex);
else
    error('Invalid option string passed into stim_*.m function');
end

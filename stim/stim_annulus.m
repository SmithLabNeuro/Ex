function stim_annulus(optstr,w,objID,arg)
%function stim_circulargrid(optstr,w,objID,arg)
%
% showex helper function for 'circulargrid' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %i %i %i');
    % arguments: (1) FrameCount
    %            (2) x position (center)
    %            (3) y position (center)
    %            (4) outer radius
    %            (5) thickness
    %            (6) color (R)
    %            (7) color (G)
    %            (8) color (B)
    stimname = mfilename;
    xPos = a(2);
    yPos = -a(3);
    
    % double check for integers and less than 0 radii
    outerRad = round(a(4));
    innerRad = round(outerRad - a(5));
    
    if outerRad <= 0
        outerRad = 0;
        innerRad = 0;
    end
    if innerRad <= 0
        innerRad = 0;
    end
    
    
    objects{objID} = struct('type', stimname(6:end),'frame',0,'fc',a(1),'x',xPos, ...
        'y',yPos,'outerRad',outerRad,'innerRad',innerRad,'col',a(6:end));
elseif strcmp(optstr,'display')
    outerPos = [sv.midScreen + [objects{objID}.x objects{objID}.y] - objects{objID}.outerRad,sv.midScreen + [objects{objID}.x objects{objID}.y] + objects{objID}.outerRad];
    Screen(w,'FillOval',objects{objID}.col,outerPos);
    
    innerPos = [sv.midScreen + [objects{objID}.x objects{objID}.y] - objects{objID}.innerRad,sv.midScreen + [objects{objID}.x objects{objID}.y] + objects{objID}.innerRad];
    Screen(w,'FillOval',sv.bgColor,innerPos);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end
function stim_arc(optstr,w,objID,arg)
%function stim_oval(optstr,w,objID,arg)
%
% showex helper function for 'arc' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %i %i %i %i %f %i %i');
    % arguments: (1) frameCount
    %            (2) x position
    %            (3) y position
    %            (4) radius
    %            (5) color, R
    %            (6) color, G
    %            (7) color, B
    %            (8) color, alpha   -13/3/19 HS for transparency
    %            (9) startAngle
    %            (10) arcAngle
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'x',a(2), ...
        'y',-a(3),'rad',a(4),'col',a(5:8),'startAng',a(9),'arcAng',a(10));
elseif strcmp(optstr,'display')
    targetPos = [sv.midScreen + [objects{objID}.x objects{objID}.y] - objects{objID}.rad,sv.midScreen + [objects{objID}.x objects{objID}.y] + objects{objID}.rad];
    Screen(w,'FillArc',objects{objID}.col,targetPos,objects{objID}.startAng,objects{objID}.arcAng);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

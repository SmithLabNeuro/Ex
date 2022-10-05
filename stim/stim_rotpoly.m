function stim_rotpoly(optstr,w,objID,arg)
%function stim_rect(optstr,w,objID,arg)
%
% showex helper function for 'rect' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %i %i %i %i %i %i %i %i %i');
    % arguments: (1) frameCount
    %            (2) x position
    %            (3) y position
    %            (4) x half-width
    %            (5) y half-height
    %            (6) angle in degrees (0 is to right)
    %            (7) color, R
    %            (8) color, G
    %            (9) color, B
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'x',a(2), ...
        'y',-a(3),'width',a(4),'height',a(5),'angle',-a(6),'col',a(7:9));
elseif strcmp(optstr,'display')
    rectCoords = [-objects{objID}.width objects{objID}.height;
             objects{objID}.width objects{objID}.height;
             objects{objID}.width -objects{objID}.height;
             -objects{objID}.width -objects{objID}.height];
    rotMat = [cosd(objects{objID}.angle) -sind(objects{objID}.angle);
            sind(objects{objID}.angle) cosd(objects{objID}.angle)];
    targetPtsInit = rectCoords * rotMat';
    targetPts = sv.midScreen + [objects{objID}.x objects{objID}.y] + targetPtsInit;
    Screen(w,'FillPoly',objects{objID}.col,targetPts);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

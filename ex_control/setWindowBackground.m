function setWindowBackground(w)
% function setWindowBackground(w)
% draws a cross hairs on an arbitrary window.  needs to be called whenever
% the voltage or eye window is redrawn with calibration dots
%
% Modified 31Oct2012 to add polar coordinate reference frame to eye window
% background -Adam C. Snyder adam@adamcsnyder.com
%
global wins params

gray=(WhiteIndex(0)+BlackIndex(0))/2;
polarColor = [170;170;170];
radii = 5:5:25; %in degrees of visual angle
radii = radii(:)'; %ensure row vector;
thetas = pi.*(0:45:359)./180; %inside the parentheses is degrees


[width height] = Screen('WindowSize',w);

if all([width height]==wins.voltageSize(3:4)-wins.voltageSize(1:2))
    windowType = 'voltage';
elseif all([width height]==wins.eyeSize(3:4)-wins.eyeSize(1:2))
    windowType = 'eye';
else
    windowType = 'other';
end;

Screen(w,'FillRect',gray);
switch windowType
    case 'voltage'
        %Note: Actually, don't show the polar coordinate frame for the
        %voltage window, because it would really have to be aligned to the
        %calibration dots to be meaningful, and there are other issues that
        %make it a bigger pain in the butt to implement than for the eye
        %window. -ACS 31Oct2012
%         pixPerDeg = [1 1]; ...deg2pix(1).*wins.pixelsPerMV;
%         xRad = radii.*pixPerDeg(1); yRad = radii.*pixPerDeg(2);
%         rects = bsxfun(@plus,[-xRad;yRad;xRad;-yRad],[wins.midV wins.midV]'); %make rectangle edges
%         Screen(w,'FrameOval',repmat(polarColor,1,numel(radii)),rects);
    case 'eye'
        pixPerDeg = deg2pix(1).*wins.pixelsPerPixel;
        xRad = radii.*pixPerDeg(1); yRad = radii.*pixPerDeg(2);
        rects = bsxfun(@plus,[-xRad;yRad;xRad;-yRad],[wins.midE wins.midE]'); %make rectangle edges
        Screen(w,'FrameOval',repmat(polarColor,1,numel(radii)),rects);
        for thx = 1:numel(thetas)
            Screen(w,'DrawLine',polarColor,...
                wins.midE(1),...
                wins.midE(2),...
                wins.midE(1)+max(xRad).*cos(thetas(thx)),...
                wins.midE(2)+max(yRad).*sin(thetas(thx)));
        end; 
        %add a circle for subject RF (-ACS): 
        %Note: at the moment only one circle is supported. In the future it
        %would be nice to have multiple circles. -ACS
        showRF =    isfield(params,'rfX') & ...
                    isfield(params,'rfY') & ...
                    isfield(params,'rfRad');
        if showRF,
            rfRad = params.rfRad;
            rfX = params.rfX; rfY = params.rfY;
            rect = [rfX-rfRad,rfY-rfRad,rfX+rfRad,rfY+rfRad].*[wins.pixelsPerPixel wins.pixelsPerPixel]+[wins.midE wins.midE];
            Screen(w,'FrameOval',polarColor,rect');
        end;
    otherwise
        %do nothing
end
Screen(w,'DrawLine',[255 255 255],0,height/2,width,height/2);
Screen(w,'DrawLine',[255 255 255],width/2,0,width/2,height);
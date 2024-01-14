function stim_radialcheck2(optstr,w,objID,arg)
%function stim_radialcheck2(optstr,w,objID,arg)
%
% showex helper function for 'radialcheck' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %i %i %i %i %i %i %i');
    % arguments: (1) FrameCount
    %            (2) dwell (frames between rotation angle changes
    %            (3) screenYpix
    %            (4) rcircles
    %            (5) tcircles
    %            (6) alpha (0-255) - not used?
    %            (7) stimX
    %            (8) stimY
    %            (9) color option: 0 for black and white only; 4 for
    %            purple/green, 1 for red/green, 2 for blue-yellow, 3 for
    %            orange-blue/green
    %            (10) rotation angle for each step of the checkerboard (starts at zero)
     [xCenter, yCenter] = RectCenter(sv.screenRect);
    screenNumber = max(Screen('Screens'));
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    grey = white / 2;
    stimname = mfilename;
    dwell = a(2); if dwell==0, dwell=1; end; % protect from divide by zero)
    screenYpix = a(3);
    rcycles = a(4);
    tcycles = a(5);
    stimX = a(7);
    stimY = a(8);
    rotAngle = a(10);

    xylim = 2 * pi * rcycles;
    [x, y] = meshgrid(-xylim: 2 * xylim / (screenYpix - 1): xylim,...
        -xylim: 2 * xylim / (screenYpix - 1): xylim);
    at = atan2(y, x);
    checks = ((1 + sign(sin(at * tcycles) + eps)...
        .* sign(sin(sqrt(x.^2 + y.^2)))) / 2) * (white - black) + black;
    circle = x.^2 + y.^2 <= xylim^2;
    checks = uint8(circle .* checks + grey * ~circle);
    
    not_checks=255-checks;
    cblank=zeros(size(checks),'uint8');
    cblank(checks==grey)=grey;
    ctrans=255*ones(size(checks),'uint8');
    ctrans(checks~=grey)=a(6);
    
    if a(9)==0
        %black-white
        color_checks{1} = cat(3,checks,checks,checks,ctrans);
    elseif a(9)==1
        %red-green
        color_checks{2} = cat(3,not_checks,checks,cblank);
    elseif a(9)==2
        %blue-yellow
        color_checks{3} = cat(3,checks,checks,not_checks);
    elseif a(9)==3
        %orange-bluegreen
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.9*255;
        color_checks{4}=cat(3,not_checks,ctemp1,checks);
    elseif a(9)==4
        %purple-green
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.5*255;
        ctemp2=checks;
        ctemp2(ctemp2==0)=0.2*255;
        color_checks{5}=cat(3,ctemp1,ctemp2,not_checks);
    end
    
    baseRect = [0 0 screenYpix screenYpix];
    dstRects(:, 1) = CenterRectOnPointd(baseRect, xCenter +stimX,yCenter+stimY);
    
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'dwell',dwell,'checks',color_checks{1},'position',dstRects);
 
    objects{objID}.radialCheckTexture = Screen('MakeTexture',w,objects{objID}.checks);
    objects{objID}.rotAngle = rotAngle; 
elseif strcmp(optstr,'display')
    rotationIdx = floor(objects{objID}.frame/objects{objID}.dwell);
    Screen('DrawTexture',w,objects{objID}.radialCheckTexture,[],objects{objID}.position,objects{objID}.rotAngle*rotationIdx);
elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.radialCheckTexture);
else
    error('Invalid option string passed into stim_*.m function');
end

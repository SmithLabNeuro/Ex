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
maxColors=5; % # of color options implemented below

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
    grey = ceil(white / 2);
    stimname = mfilename;
    dwell = a(2); if dwell==0, dwell=1; end; % protect from divide by zero)
    screenYpix = a(3);
    rcycles = a(4);
    tcycles = a(5);
    stimX = a(7);
    stimY = a(8);
    chkCol = a(9);
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
    not_checks(not_checks==127)=grey; % fix the gray
    cblank=uint8(circle)-1;
    cblank(cblank==1)=grey;
    ctrans=255*ones(size(checks),'uint8');
    ctrans(checks~=grey)=a(6);
    
    if chkCol==0
        %black-white
        color_checks{1} = cat(3,checks,checks,checks,ctrans);
    elseif chkCol==1
        %red-green
        color_checks{2} = cat(3,not_checks,checks,cblank,ctrans);
    elseif chkCol==2
        %blue-yellow
        color_checks{3} = cat(3,checks,checks,not_checks,ctrans);
    elseif chkCol==3
        %orange-bluegreen
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.9*255;
        color_checks{4}=cat(3,not_checks,ctemp1,checks,ctrans);
    elseif chkCol==4
        %purple-green
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.5*255;
        ctemp2=checks;
        ctemp2(ctemp2==0)=0.2*255;
        color_checks{5}=cat(3,ctemp1,ctemp2,not_checks,ctrans);
    elseif chkCol==-1 % make all the checks
        color_checks{1} = cat(3,checks,checks,checks,ctrans);
        color_checks{2} = cat(3,not_checks,checks,cblank,ctrans);
        color_checks{3} = cat(3,checks,checks,not_checks,ctrans);
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.9*255;
        color_checks{4}=cat(3,not_checks,ctemp1,checks,ctrans);
        ctemp1=checks;
        ctemp1(ctemp1==0)=0.6*255;
        ctemp1(ctemp1==255)=0.5*255;
        ctemp2=checks;
        ctemp2(ctemp2==0)=0.2*255;
        color_checks{5}=cat(3,ctemp1,ctemp2,not_checks,ctrans);
    else
        error('Invalid checkerboard color parameter');
    end
    
    baseRect = [0 0 screenYpix screenYpix];
    dstRects(:, 1) = CenterRectOnPointd(baseRect, xCenter +stimX,yCenter+stimY);
    
    if chkCol>=0
        objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'dwell',dwell,'checks',color_checks{chkCol+1},'position',dstRects);
        objects{objID}.radialCheckTexture = struct();
        objects{objID}.rotAngle = rotAngle;
        objects{objID}.chkCol=chkCol;
        objects{objID}.radialCheckTexture(1).tex = Screen('MakeTexture',w,objects{objID}.checks);
    elseif chkCol==-1
        objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'dwell',dwell,'position',dstRects);
        objects{objID}.radialCheckTexture = struct();
        objects{objID}.rotAngle = rotAngle;
        objects{objID}.chkCol=chkCol;          
        objects{objID}.checks = color_checks;
        for I=1:maxColors
            objects{objID}.radialCheckTexture(I).tex = Screen('MakeTexture',w,color_checks{I});
        end
    end
elseif strcmp(optstr,'display')
  rotationIdx = floor(objects{objID}.frame/objects{objID}.dwell);
  if objects{objID}.chkCol>=0
    Screen('DrawTexture',w,objects{objID}.radialCheckTexture(1).tex,[],objects{objID}.position,objects{objID}.rotAngle*rotationIdx);
  elseif objects{objID}.chkCol==-1
    colorIdx = mod(rotationIdx+1,maxColors);
    if colorIdx == 0
      colorIdx = maxColors;
    end
    Screen('DrawTexture',w,objects{objID}.radialCheckTexture(colorIdx).tex,[],objects{objID}.position,objects{objID}.rotAngle*rotationIdx);
  end
elseif strcmp(optstr,'cleanup')
  for I=1:numel(objects{objID}.radialCheckTexture)
    Screen('Close',objects{objID}.radialCheckTexture(I).tex);    
  end
  %
  %if objects{objID}.chkCol>=0
  %      Screen('Close',objects{objID}.radialCheckTexture(1).tex);
  %  elseif objects{objID}.chkCol==-1
  %      for I=1:maxColors
  %          Screen('Close',objects{objID}.radialCheckTexture(I).tex);
  %      end
  %  end
else
    error('Invalid option string passed into stim_*.m function');
end

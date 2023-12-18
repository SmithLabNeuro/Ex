function stim_circulargrid(optstr,w,objID,arg)
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
    a = sscanf(arg,'%i %i %i %i %f %i %i %i %i');
    % arguments: (1) FrameCount
    %            (2) screenYpix
    %            (3) rcircles
    %            (4) tcircles
    %            (5) alpha (0-255)
    %            (6) stimX
    %            (7) stimY
    %            (8) color option: 0 for black and white only; 1 for random
    %            color
    %            (9) cue type (1,2,3 => S, M, L)
    [xCenter, yCenter] = RectCenter(sv.screenRect);
    screenNumber = max(Screen('Screens'));
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    grey = white / 2;
    stimname = mfilename;
    screenYpix = a(2);
    rcycles = a(3);
    tcycles = a(4);
    stimX = a(6);
    stimY = -a(7);
    xylim = 5;
  
    cuetype = a(9);
    if(cuetype == 1)
        checks = [0,1,1,1,1;...
                  0,1,0,0,0;...
                  0,1,1,1,0;...
                  0,0,0,1,0;...
                  1,1,1,1,0];
    elseif(cuetype == 2)
        checks = [1,0,0,0,1;...
                  1,1,0,1,1;...
                  1,0,1,0,1;...
                  1,0,0,0,1;...
                  1,0,0,0,1];
    elseif(cuetype == 3)
        checks = [1,0,0,0,0;...
                  1,0,0,0,0;...
                  1,0,0,0,0;...
                  1,1,1,1,1;...
                  1,1,1,1,1];
    elseif(cuetype == 4)
        checks = [1,1,1,1,1;...
                  0,0,1,0,0;...
                  0,0,1,0,0;...
                  1,1,1,0,0;...
                  1,1,1,0,0];
    end
    checks = 1-checks;
    checks = checks.*255;
    
    not_checks=255-checks;
    cblank=zeros(size(checks));
    cblank(checks==grey)=grey;
    ctrans=255*ones(size(checks));
    ctrans(checks~=grey)=a(5);
    
    %black-white
    color_checks{1} = cat(3,checks,checks,checks);

    %red-green
    color_checks{2} = cat(3,not_checks,checks,cblank);
    
    %blue-yellow
    color_checks{3} = cat(3,checks,checks,not_checks);
    
    %orange-bluegreen
    ctemp1=checks;
    ctemp1(ctemp1==0)=0.6*255;
    ctemp1(ctemp1==255)=0.9*255;
    color_checks{4}=cat(3,not_checks,ctemp1,checks);
    
    %purple-green
    ctemp1=checks;
    ctemp1(ctemp1==0)=0.6*255;
    ctemp1(ctemp1==255)=0.5*255;
    ctemp2=checks;
    ctemp2(ctemp2==0)=0.2*255;
    color_checks{5}=cat(3,ctemp1,ctemp2,not_checks);
    
    
    for i = 1:length(color_checks)
        for n = 1:size(checks,1)
            for m = 1:size(checks,2)
                if checks(n,m)==grey
                    color_checks{i}(n,m,1:3) = sv.bgColor;
                end
            end
        end
        color_checks{i} = cat(3,color_checks{i},ctrans);
    end
    baseRect = [0 0 screenYpix screenYpix];
    dstRects(:, 1) = CenterRectOnPointd(baseRect, xCenter +stimX,yCenter+stimY);
    % Create a separate oval rectangle that's larger than the rectangle for
    % the grid
    baseOvalRect = [0, 0, screenYpix*sqrt(2), screenYpix*sqrt(2)];
    ovalRects(:,1) = CenterRectOnPointd(baseOvalRect, xCenter+stimX,yCenter+stimY);
    if numel(a)<8 || a(8) == 0
        objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), 'col', white, 'checks',color_checks{1},'position',dstRects, 'ovalPosition', ovalRects);
    else
        objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), 'col', white, 'checks',color_checks{randi(5)},'position',dstRects, 'ovalPosition', ovalRects);
    end
elseif strcmp(optstr,'display')
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    radialCheckerboardTexture  = Screen('MakeTexture', w, objects{objID}.checks);
    Screen('FillOval', w, objects{objID}.col, objects{objID}.ovalPosition);
    Screen('DrawTexture', w, radialCheckerboardTexture,[],objects{objID}.position);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end
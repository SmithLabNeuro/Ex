function stim_squarecheck(optstr,w,objID,arg)
%function stim_radialcheck(optstr,w,objID,arg)
%
% showex helper function for 'oval' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %f %f %f %f %f %f %f %f %i %i %i %i %i');
    % arguments: (1) FrameCount
    %            (2) rgrid
    %            (3) ngrid
    %            (4) color 1 R
    %            (5) color 1 G
    %            (6) color 1 B
    %            (7) color 1, alpha
    %            (8) color 2 R
    %            (9) color 2 G
    %            (10) color 2 B
    %            (11) color 2, alpha
    %            (12) seed
    %            (13) refresh
    %            (14) stimX
    %            (15) stimY
    %            (16) aperture
    [xCenter, yCenter] = RectCenter(sv.screenRect);
    stimname = mfilename;
    rgrid = a(2);
    ngrid = a(3);
    stim_x = a(14);
    stim_y = a(15);
    if numel(a)<12
        seed = 0;
    else
        seed = a(12);
    end
    if numel(a)<13
        refresh = 0;
    else
        refresh = a(13);
    end
    if numel(a)<16
        aperture = -1;
    else
        aperture = a(16);
    end
    if numel(a)>=11
        color1 = a(4:7);
        color2 = a(8:11);
    elseif numel(a) >=7
        color1 = a(4:7);
        color2 = [0,0,0,0.5];
    else
        color1 = [255,255,255,0.5];
        color2 = [0,0,0,0.5];
    end
    baseRect = [0 0 rgrid rgrid];
    if mod(ngrid,2)==1
        [xPos, yPos] = meshgrid(-(ngrid-1)/2:1:(ngrid-1)/2, -(ngrid-1)/2:1:(ngrid-1)/2);
        ch=repmat(eye(2),(ngrid+1)/2,(ngrid+1)/2);
        ch = ch(1:end-1,1:end-1);
    else
        [xPos, yPos] = meshgrid(-ngrid/2:1:ngrid/2-1, -ngrid/2:1:ngrid/2-1);
        ch=repmat(eye(2),ngrid/2,ngrid/2);
    end
    [s1, s2] = size(xPos);
    numSquares = s1 * s2;
    
    xPos = reshape(xPos, 1, numSquares).* rgrid;
    yPos = reshape(yPos, 1, numSquares).* rgrid;
    ch1= reshape(ch,size(ch,1)*size(ch,2),1);
    if aperture ~= -1
        % check the center position of each square. compare it to the
        % aperture size in pixels, if it's greater, then delete that square
        removeList = [];
        for i = 1:length(xPos)
            if xPos(i)^2+yPos(i)^2>aperture^2
                removeList = [removeList,i];
            end
        end
    end
    xPos = xPos + xCenter + stim_x;
    yPos = yPos + yCenter - stim_y;

    for i=1:size(ch1,1)
        if ch1(i)
            colors(:,i)=color1';
        else
            colors(:,i)=color2';
        end
    end
    r = RandStream.create('mrg32k3a','seed',seed);
    idx = randperm(r,length(ch1));
    if seed>0
        colors = colors(:,idx);
    end
    allRects = nan(4, 3);
    for i = 1:numel(xPos)
        allRects(:, i) = CenterRectOnPointd(baseRect,...
            xPos(i), yPos(i));
    end
    if aperture ~= -1
        allRects(:,removeList)=[];
        colors(:,removeList)=[];
    end
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1),'allRects',allRects,'colors',colors,'seed',seed,'refresh',refresh);
elseif strcmp(optstr,'display')
    
    if objects{objID}.refresh == 1
        rmseed = randi(10);
        r = RandStream.create('mrg32k3a','seed',rmseed);
        idx = randperm(r,size(objects{objID}.colors,2));
        objects{objID}.colors = objects{objID}.colors(:,idx);
    end
    % updates the colors randperm
    % uses 'r' seed to make these new patterns.
    
    Screen('FillRect',w,objects{objID}.colors,objects{objID}.allRects);
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

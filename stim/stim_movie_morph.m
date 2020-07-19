function stim_movie_morph2(optstr,w,objID,arg)
%function stim_movie_morph(optstr,w,objID,arg)
%
% showex helper function for 'movie_morph' stim class that can
% alter the movie frames when the movie is read. This allows you to
% have a base movie (or image) and then alter it on the fly at read time.
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    [fileName argsRest] = strtok(arg);
    a = sscanf(argsRest,'%i %i %i %i %i %i %i');
    
    % arguments: (1) frameCount
    %            (2) dwell
    %            (3) start frame
    %            (4) x position
    %            (5) y position
    %            (6) aperture (0 for original dimensions, specify
    %            pixel diameter if you want a circle)
    %            (7) morph type
    %            (8) morph angle
    %            (9) width  % optional rescaling values
    %            (10) height % optional rescaling values
    
    dwell = a(2);
    startFrame = a(3);
    xCenter = a(4);
    yCenter = a(5) * -1; % flip y coordinate so '-' is down
    apDim = a(6);
    morphType = a(7);
    morphAngle = a(8);
    
    if numel(a) > 9
        rescaleWidth = a(9);
        rescaleHeight = a(10);
    elseif numel(a) == 9
        rescaleWidth = a(9);
        rescaleHeight = rescaleWidth;
    end
    
    vars = load([sv.localDir,filesep,fileName]); %added sv.localDir, now prefix/suffix combos will be relative to that location. This directory can be set by rig/subject XMLs. -ACS 13Sep2013
    % PLS 07/05/19 rename "img" variable to "mov"
    if ~isfield(vars,'mov')
        fields = struct2cell(vars);
        vars.mov{1} = uint8(fields{1}.*255);
    end
    
    % rescale if requested
    if numel(a) > 8
        for I=1:numel(vars.mov)
            vars.mov{I}=imresize(vars.mov{I},[rescaleWidth rescaleHeight],'nearest');
        end
    end
    
    movieSize = size(vars.mov{1});
    
    % Mask with a circular aperture
    if apDim == 0
        apDim = -1;
    end
    
    % movieXSize = movieSize(1);
    %movieYSize = movieSize(2);
    %movieXCtr = round(movieXSize/2);
    %movieYCtr = round(movieYSize/2);
    %for F=1:length(vars.mov)
    %    disp(['masking frame ',num2str(F)]);
    %    for I=1:movieXSize
    %        for J=1:movieYSize
    %            x = I-movieXCtr;
    %            y = J-movieYCtr;
    %            if (x*x + y*y > apDim)
    %                vars.mov{F}(I,J) = 0; % use bg value instead
    %            end
    %        end
    %    end
    %end
    
    maskDimX = (movieSize(1)-1)/2;
    maskDimY = (movieSize(2)-1)/2;
    
    % Create a single gaussian transparency mask and store it to a texture:
    mask=ones(movieSize(1),movieSize(2), 2) * mean(sv.bgColor);
    
    [x,y]=meshgrid(-1*maskDimX:1*maskDimX,-1*maskDimY:1*maskDimY);
    
    mask(:, :, 2)=sv.white * (sqrt(x.^2+y.^2)' > apDim);
    
    maskTex=Screen('MakeTexture', w, mask);
    
    % Now morph the images however you want
    if morphType > 0
        for F=1:length(vars.mov)
            if morphType == 1 % L*a*b* color shift
                cform = makecform('srgb2lab');
                invcform = makecform('lab2srgb');
                labImg = applycform(im2double(vars.mov{F}),cform);
                theta = deg2rad(morphAngle); 
                rotmat = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                colors = reshape(labImg(:,:,2:3),[],2);
                newColors = colors*rotmat;
                newColors = reshape(newColors,movieSize(1),movieSize(2),2);
                labAnew = cat(3,labImg(:,:,1),newColors);
                vars.mov{F} = uint8(applycform(labAnew,invcform).*255); %new image
            elseif morphType == 2 % orientation shift
                vars.mov{F} = imrotate(vars.mov{F},morphAngle,'bilinear','crop');
            elseif morphType == 3 % zoom
                origArea = movieSize(1)*movieSize(2);
                origX = movieSize(1);
                origY = movieSize(2);
                center = [round(origX/2) round(origY/2)];
                newArea = origArea*(1-morphAngle/100);
                newX = sqrt(newArea*origX/origY);
                newY = newArea/newX;
                vars.mov{F} = vars.mov{F}(center(1)-floor(newY/2):center(1)+ceil(newY/2),center(2)-floor(newX/2):center(2)+ceil(newX/2),:);
                vars.mov{F}=imresize(vars.mov{F},[movieSize(1) movieSize(2)],'nearest');
            elseif morphType == 4 %L*a*b (L)Lightness shift
                cform = makecform('srgb2lab');
                invcform = makecform('lab2srgb');
                labImg = applycform(im2double(vars.mov{F}),cform);
                labImg(:,:,1) = labImg(:,:,1)+morphAngle;
                labImg(labImg(:,:,1)>100) = 100;
                labImg(labImg(:,:,1)<0) = 0;                
                vars.mov{F} = uint8(applycform(labImg,invcform).*255); %new image
            elseif morphType == 5 %hsv (h)hue shift
                hsvimg = rgb2hsv(vars.mov{F});
                hsvimg(:,:,1) = hsvimg(:,:,1) + morphAngle*(1/360); %increase hue
                hsvimg(hsvimg>1) = hsvimg(hsvimg>1)-1;
                hsvimg(hsvimg<0) = abs(hsvimg(hsvimg<0));
                vars.mov{F} = hsv2rgb(hsvimg).*255;
            elseif morphType == 6 % hsv (s)saturation
                hsvimg = rgb2hsv(vars.mov{F});
                hsvimg(:,:,2) = hsvimg(:,:,2) + morphAngle/100; %increase saturation
                hsvimg(hsvimg>1) = 1;
                hsvimg(hsvimg<0) = 0;
                vars.mov{F} = hsv2rgb(hsvimg).*255;
            elseif morphType == 7 % hsv (v)value shift
                hsvimg = rgb2hsv(vars.mov{F});
                hsvimg(:,:,3) = hsvimg(:,:,3) + morphAngle/100; %increase value
                hsvimg(hsvimg>1) = 1;
                hsvimg(hsvimg<0) = 0;
                vars.mov{F} = hsv2rgb(hsvimg).*255;
            end
        end
    end
    
    % MAS 06/27/13 - removed the "-1" from
    % dstRect/srcRect. Also made the code use the first
    % two values in movieSize for when a 3D (color)
    % movie is in use
    dstRect=[0 0 movieSize(2:-1:1)];
    dstRect=CenterRect(dstRect, sv.screenRect) + [xCenter yCenter xCenter yCenter];
    srcRect=[0 0 movieSize(2:-1:1)];

    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1), ...
        'startFrame',startFrame, 'dwell',dwell, ...
        'srcRect',srcRect, 'dstRect',dstRect,'mask',maskTex);
    
    objects{objID}.mov = vars.mov;
elseif strcmp(optstr,'display')
    frameNum = floor(objects{objID}.frame/objects{objID}.dwell) + objects{objID}.startFrame;
    
    if frameNum > length(objects{objID}.mov)
        frameNum = mod(frameNum-1,length(objects{objID}.mov))+1;
    end
    
    tex=Screen('MakeTexture', w, objects{objID}.mov{frameNum});
    Screen('DrawTexture', w, tex,objects{objID}.srcRect,objects{objID}.dstRect);
    Screen('DrawTexture',w,objects{objID}.mask,objects{objID}.srcRect,objects{objID}.dstRect);
    Screen('Close',tex);
elseif strcmp(optstr,'cleanup')
    Screen('Close',objects{objID}.mask);    
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

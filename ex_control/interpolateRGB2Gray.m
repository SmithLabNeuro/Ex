function rgbvec = interpolateRGB2Gray(startRGB,grayscalar,mixamount,hsvflag)
%if ~exist('hsvflag','var');hsvflag = 0;end

if hsvflag == 1
    %% convert rgb endpoints to hsv
    normstartRGB = startRGB/255;
    normgrayscalar = [1 1 1]*grayscalar/255;
    starthsv = rgb2hsv(normstartRGB);
    endhsv = rgb2hsv(normgrayscalar);

    %% interpolate in hsv space
    hsvgrad = zeros(length(mixamount),3);
    hsvgrad(:,1) = repmat(starthsv(1),length(mixamount),1);
    for n = 2:3
        hsvgrad(:,n) =(endhsv(n)-starthsv(n))*mixamount + starthsv(n);
    end

    %% convert to rgb
    rgbvec = round(hsv2rgb(hsvgrad)*255);
else
    if length(grayscalar) == 1
        endrgb = [1 1 1]*grayscalar;
    else
        endrgb = grayscalar;
    end
    rgbvecpre = zeros(length(mixamount),3);
    for n = 1:3
        rgbvecpre(:,n) =(endrgb(n)-startRGB(n))*mixamount + startRGB(n);
    end
    rgbvec = round(rgbvecpre);
end
end
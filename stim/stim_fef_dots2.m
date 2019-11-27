function stim_fef_dots2(optstr,w,objID,arg)
%function stim_fef_dots(optstr,w,objID,arg)
%
% showex helper function for 'fef_dots' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i %i %i %f %i %i %i %i %i %i %i %i');
    
    % arguments: (1) frameCount
    %            (2) random seed
    %            (3) number of dots per frame
    %            (4) dot size
    %            (5) dwell
    %            (6) center x
    %            (7) center y
    %            (8) x radius
    %            (9) y radius
    %            (10) color, R
    %            (11) color, G
    %            (12) color, B
    
    frameCount = a(1);
    seed = a(2);
    numDots = a(3);
    dotRad = a(4)*3;
    dwell = a(5);
    
    xMin = (a(6) - a(8));
    xMax = (a(6) + a(8));
    yMin = (a(7) - a(9));
    yMax = (a(7) + a(9));
    
    dotPositions = zeros(4,numDots,ceil(frameCount/dwell));
    
    r = RandStream.create('mrg32k3a','seed',seed);
    
    for i = 1:size(dotPositions,3)
        if numDots==1
            xPos = randi(r,xMax-xMin+1,numDots,1)+xMin-1;
            yPos = randi(r,yMax-yMin+1,numDots,1)+yMin-1;
        else
            
            if(mod(i,2)==0)
                xPos = [randi(r,floor(a(8)/2),floor(numDots/2),1)+xMin-1; randi(r,a(8)/2,ceil(numDots/2),1)+a(8)/2*2+xMin-1];
            else
                xPos = [randi(r,floor(a(8)/2),floor(numDots/2),1)+a(8)/2+xMin-1; randi(r,a(8)/2,ceil(numDots/2),1)+a(8)/2*3+xMin-1];
            end
        end
        
        yPos = randi(r,yMax-yMin+1,numDots,1)+yMin-1;
        
        %%%%%% begin scalediam function %%%%%%
        % xPos, yPos and dotRad were in
        % units of degrees - is this necessary?
        %E = sqrt((xPos./ppd).^2 + (yPos./ppd).^2);
        %X = log10(E) - 1.5;
        %equality = 0.8124 + (0.5324 .* X) + (0.0648 * X.^2) + (0.0788 * X.^3);
        %N = 10.^equality;
        %M = 1.0 ./ N;
        %alldotRad = sqrt(((dotRad./ppd)^2)./M); % dot radii
        %alldotRad = alldotRad .* ppd;
        %%%%%% end scalediam function %%%%%%
        % this catches the case where the dot is at 0,0
        %if (sum(isnan(alldotRad)) > 0)
        %    alldotRad(find(isnan(alldotRad))) = 0;
        %    display(i)
        %end
        
        dotPositions(1,:,i) = xPos - dotRad;
        dotPositions(2,:,i) = yPos - dotRad;
        dotPositions(3,:,i) = xPos + dotRad;
        dotPositions(4,:,i) = yPos + dotRad;
    end
    
    dotPositions([1 3],:,:) = dotPositions([1 3],:,:) + sv.midScreen(1);
    dotPositions([2 4],:,:) = dotPositions([2 4],:,:) + sv.midScreen(2);
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',frameCount,'dp',dotPositions, ...
                            'color',a(10:12),'dwell',a(5));
elseif strcmp(optstr,'display')
    Screen(w,'FillRect',objects{objID}.color,objects{objID}.dp(:,:,1+floor(objects{objID}.frame/objects{objID}.dwell)));
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end


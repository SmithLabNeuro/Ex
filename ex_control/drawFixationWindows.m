function drawFixationWindows(fixX,fixY,rAll,varargin)
% function drawFixationWindows(fixX,fixY,r)
% function drawFixationWindows(fixX,fixY,r,winColors)
% 
% ex helper function. draws the fixation windows to the screen. used by
% functions like waitForFixation, waitForMS, and waitForSlave, since they
% require the fixation window parameters as input.
%
% If a fourth argument, winColors, should be a N x 3 matrix where N is the
% number of windows (and the 3 values are the RGB window colors from
% 0-255). If only one RGB value is given, it applies to all windows. If
% winColors argument is not provided, the windows are yellow by default.
%
% fixX, fixY and r can be vectors for multiple windows
%
%Modified 28Mar2012 by Adam Snyder to support multiple windows (and colors)

global calibration wins;
   
    Screen('CopyWindow',wins.voltageBG,wins.voltage,wins.voltageDim,wins.voltageDim);
    Screen('CopyWindow',wins.eyeBG,wins.eye,wins.eyeDim,wins.eyeDim);
    if nargin > 0
        fixY = -fixY; %flip y coordinate because PTB's native coordinate space has negative up, but Ex uses negative down. -ACS13Mar2012
        numWindows = unique([length(fixX) length(fixY) size(rAll,2)]); 
        assert(numel(numWindows)==1,'Fixation window parameters X, Y and R must be same size');
        if nargin > 3
            winColors = varargin{1};
            if size(winColors,1)>1
                assert(size(winColors,1)==numWindows,'If more than one color is provided then number of colors must equal the number of windows');
            else
                winColors = repmat(winColors(:)',numWindows,1);
            end;
        else
            winColors = repmat([255 255 0],numWindows,1);
        end;                
        for i = 1:numWindows
            r = rAll(~isnan(rAll(:, i)), i);
            if size(r,1)==1 %indicates a circle is desired
%                 fixationWindow = r.*[-1 -1; -1 1; 1 1; 1 -1] + repmat([fixX(i) fixY(i)],4,1); 
                fixationWindow = r.*[-1 -1 1 1]' + [fixX(i) fixY(i) fixX(i) fixY(i)]'; 
                Screen('FrameOval',wins.eye,winColors(i,:),fixationWindow.*repmat(wins.pixelsPerPixel,1,2)'+repmat(wins.midE,1,2)');
                fixationWindow = reshape(fixationWindow,2,[])'; %put X values in one column and Y values in another...
                vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
                vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
                vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
                vPoints = reshape(vPoints',[],1); %put back to column form
%                 disp(vPoints);
                Screen('FrameOval',wins.voltage,winColors(i,:),vPoints);
            elseif size(r,1)==2
                fixationWindow = bsxfun(@times,r',[-1 -1; -1 1; 1 1; 1 -1]) + repmat([fixX(i) fixY(i)],4,1); 
                Screen('FramePoly',wins.eye,winColors(i,:),fixationWindow.*repmat(wins.pixelsPerPixel,4,1)+repmat(wins.midE,4,1));
                vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
                vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
                vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
                Screen('FramePoly',wins.voltage,winColors(i,:),vPoints);
            elseif size(r,1)==3
                % wedge
                if false%i == 1 || i == 2
                    fixationWindow = r(1).*[-1 -1 1 1]' + [fixX(i) fixY(i) fixX(i) fixY(i)]'; 
                    Screen('FrameOval',wins.eye,winColors(i,:),fixationWindow.*repmat(wins.pixelsPerPixel,1,2)'+repmat(wins.midE,1,2)');
                    fixationWindow = reshape(fixationWindow,2,[])'; %put X values in one column and Y values in another...
                    vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
                    vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
                    vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
                    vPoints = reshape(vPoints',[],1); %put back to column form
%                 disp(vPoints);
                    Screen('FrameOval',wins.voltage,winColors(i,:),vPoints); clear vPoints
                else
                    inrad = r(1); % in pixels
                    width = r(2); % width from middle of wedge to one side of wedge in degrees
                    ang = r(3); % angle of wedge middle in degrees
                    [px,py] = pol2cart(deg2rad([ang-width;ang-width;ang;ang+width;ang+width])',[1000;inrad;inrad;inrad;1000]');
                    fixationWindow = [px',-py'] ;%+ repmat([fixX(i) fixY(i)],length(px),1); 
                    Screen('FramePoly',wins.eye,winColors(i,:),fixationWindow.*repmat(wins.pixelsPerPixel,length(px),1)+repmat(wins.midE,length(px),1));
                    vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
                    vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
                    vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
                    Screen('FramePoly',wins.voltage,winColors(i,:),vPoints);
                end
                % if using an ellipse
                %                 fixationWindow = r(i).*[-1 -1; -1 1; 1 1; 1 -1] + repmat([fixX(i) fixY(i)],4,1);
%                 rotatecoord = ([cos(-1*r(3,i)) sin(-1*r(3,i));-sin(-1*r(3,i)) cos(-1*r(3,i))]);
%                 fixationWindow = round(bsxfun(@times,r(1:2,i)',[-1 -1; -1 1; 1 1; 1 -1])*rotatecoord) + repmat([fixX(i) fixY(i)],4,1);
%                 %fixationWindow = bsxfun(@times,r(:,i)',[-1 -1; -1 1; 1 1; 1 -1]) + repmat([fixX(i) fixY(i)],4,1); 
%                 Screen('FramePoly',wins.eye,winColors(i,:),fixationWindow.*repmat(wins.pixelsPerPixel,4,1)+repmat(wins.midE,4,1));
%                 vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
%                 vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
%                 vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
%                 Screen('FramePoly',wins.voltage,winColors(i,:),vPoints);
                
                
                %                 fixationWindow = [-r(1,i) -r(2,i) r(1,i) r(2,i)]' + [fixX(i) fixY(i) fixX(i) fixY(i)]'; 
%                 pixelwindow = fixationWindow.*repmat(wins.pixelsPerPixel,1,2)'+repmat(wins.midE,1,2)';
%                 Screen('glPushMatrix', wins.eye)
%                 Screen('glTranslate', wins.eye, (pixelwindow(1)+pixelwindow(3)/2), (pixelwindow(2)+pixelwindow(4)/2))
%                 Screen('glRotate', wins.eye, 45, 0, 0);
%                 Screen('glTranslate', wins.eye, -(pixelwindow(1)+pixelwindow(3)/2), -(pixelwindow(2)+pixelwindow(4)/2))
%                 Screen('FrameOval',wins.eye,winColors(i,:),pixelwindow);
%                 Screen('glPopMatrix', wins.eye)
%                 fixationWindow = reshape(fixationWindow,2,[])'; %put X values in one column and Y values in another...
%                 vPoints(:,1) = [fixationWindow ones(size(fixationWindow,1),1)] * calibration{5};
%                 vPoints(:,2) = [bsxfun(@times,fixationWindow,[1 -1]) ones(size(fixationWindow,1),1)] * calibration{6};
%                 vPoints = sort(vPoints,'ascend'); %forces corners defined to be bottom-left and upper-right -ACS 24OCT2013
%                 vPoints = reshape(vPoints',[],1); %put back to column form
% %                 disp(vPoints);
%                 Screen('FrameOval',wins.voltage,winColors(i,:),vPoints); clear vPoints
            else
                error('''r'' is expected to be either 1xNwindows (circles) or 2xNwindows (rectangles)');
            end;
            clear vPoints
        end;
    end
end
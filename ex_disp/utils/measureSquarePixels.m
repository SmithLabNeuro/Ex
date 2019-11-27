%function measureSquarePixels(sqsize)
%
% This function will display a square of a known size so that the screen
% can be adjusted to make the pixels square (important for CRTs). "sqsize"
% is the size of the square to be displayed.
%

function measureSquarePixels(sqsize)

if (nargin < 1)
    sqsize = 200;
end

% Don't need gamma correction for this, but you could use it if desired
%Screen('LoadNormalizedGammaTable',0);

% open a new window, fill with gray background
w = Screen('OpenWindow',0);
Screen(w,'FillRect',[128 128 128]);
Screen('Flip',w);

% get the screen resolution
res = Screen('Resolution',0);

% figure out the positon of the square and draw it
ctrx = res.width/2;
ctry = res.height/2;
sqoffset = sqsize/2;
sqpos = [ctrx-sqoffset ctry-sqoffset ctrx+sqoffset ctry+sqoffset];
Screen(w,'FillRect',[255 255 255],sqpos);
Screen('DrawText',w,['Square is ',num2str(sqsize), ' pixels on each side'],50,75);
Screen('DrawText',w,['Grid is 100 pixels square'],50,100);

% draw grid lines on the whole screen with 100 pixel spacing
ls = 100;
nvertlines = floor(res.width / ls);
nhorizlines = floor(res.height / ls);
hlines = [fliplr([0:-ls:-ctry]),[0+ls:ls:ctry]]+ctry
vlines = [fliplr([0:-ls:-ctrx]),[0+ls:ls:ctrx]]+ctrx
for I=1:nvertlines
    Screen('DrawLine',w,[0 0 0],vlines(I),0,vlines(I),res.height);
end
for I=1:nhorizlines
    Screen('DrawLine',w,[0 0 0],0,hlines(I),res.width,hlines(I));
end

% draw some cross-hairs
ll = 10; % length of crosshair line (/2)
tctrx = ctrx; tctry = ctry;
Screen('DrawLine',w,[255 0 0],tctrx,tctry-ll,tctrx,tctry+ll,2);
Screen('DrawLine',w,[255 0 0],tctrx-ll,tctry,tctrx+ll,tctry,2);
tctrx = ctrx-sqoffset; tctry = ctry-sqoffset; 
Screen('DrawLine',w,[255 0 0],tctrx,tctry-ll,tctrx,tctry+ll,2);
Screen('DrawLine',w,[255 0 0],tctrx-ll,tctry,tctrx+ll,tctry,2);
tctrx = ctrx+sqoffset; tctry = ctry+sqoffset; 
Screen('DrawLine',w,[255 0 0],tctrx,tctry-ll,tctrx,tctry+ll,2);
Screen('DrawLine',w,[255 0 0],tctrx-ll,tctry,tctrx+ll,tctry,2);
tctrx = ctrx+sqoffset; tctry = ctry-sqoffset; 
Screen('DrawLine',w,[255 0 0],tctrx,tctry-ll,tctrx,tctry+ll,2);
Screen('DrawLine',w,[255 0 0],tctrx-ll,tctry,tctrx+ll,tctry,2);
tctrx = ctrx-sqoffset; tctry = ctry+sqoffset; 
Screen('DrawLine',w,[255 0 0],tctrx,tctry-ll,tctrx,tctry+ll,2);
Screen('DrawLine',w,[255 0 0],tctrx-ll,tctry,tctrx+ll,tctry,2);

% refresh the screen
Screen('Flip',w);

pause; % press a key to exit
    
sca

%function testSquare(vals)
%
% This function first loads an identity (uncorrected) gamma table.
% Then it brings up a gray square that fills the screen
% and shows it, using the list of grayscale values provided (0-255)
% and waiting for a space bar press between changing values.
% Useful for testing calibration quickly.
%

function testSquare(vals,size)

w = Screen('OpenWindow',0);
LoadIdentityClut(w);
Screen(w,'FillRect',0);

Screen('Flip',w);

rgbvals = [vals;vals;vals];
for i = 1:length(vals)
    if (nargin < 2)
        Screen(w,'FillRect',rgbvals(:,i));
    else
        res = Screen('Resolution',w);
        ctrx = round(res.width/2);
        ctry = round(res.height/2);
        % draw background as mean gray
        Screen(w,'FillRect',[128 128 128]);
        % draw square on top
        Screen(w,'FillRect',rgbvals(:,i),[ctrx-size ctry-size ctrx+size ctry+size]);
    end
    
    Screen('DrawText',w,num2str(rgbvals(:,i)'),100,100);
    Screen('Flip',w);
       
    pause;
    
end


sca

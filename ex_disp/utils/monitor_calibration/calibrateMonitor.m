%function calibVals = calibrateMonitor(vals,color,size,portstring)
% OR:
%function calibVals = calibrateMonitor(vals)
%function calibVals = calibrateMonitor(vals,color)
%function calibVals = calibrateMonitor(vals,color,size)
%
% This function will display different gray levels on the screen and
% communicate with a UDT 471 photometer over the serial port to read
% luminance values. It returns the list of luminance values (in cd/m^2)
% associated with the grayscale values that were passed in.
%
% The input "vals" should be a list of luminance values from 0 to 255 to be
% displayed.
%
% The input "color" can be 'gray', 'red', 'green', or 'blue' (default is
% gray if no parameter is passed), and this determines the color to be
% calibrated. If 'other' is used, you specify your own 3xN list of rgb
% values. If no color is specified, 'gray' is the default.
%
% Size can be -1 if you want it full screen, otherwise it's the
% half-diameter of the rectangle. If no size is specified, full screen is
% the default.
%
% The serial port string for opening (e.g., '/dev/ttyS0'). If no string is
% specified, it assumes '/dev/cu.KeySerial1' on the mac and '/dev/ttyS0/'
% on Linux.
%
% e.g., lum = calibrateMonitor(0:255,'gray');
%
% NOTE: If the program hangs there's a good chance that your photometer is
% not powered on, which will cause the serial port calls to hang.
%
% NOTE: This program waits 5 seconds between luminance changes to let the
% value on the UDT 471 stabilize. You can read faster than that for
% testing, but I've found 5 seconds to be a good choice for the final
% calibration measurements.
%
% Matthew A. Smith
% Revised: 20120413

function calibVals = calibrateMonitor(vals,color,size,portstring)

if (nargin < 2)
    color = 'gray';
end

delete(instrfind);

ctype = computer('arch');

if (nargin < 4)
    if (strcmp(ctype,'maci') || strcmp(ctype,'maci64'))
        % For Mac Pro with Keyspan USB Serial Installed
        display('Mac detected - opening serial port /dev/cu.KeySerial1');
        out = serial('/dev/cu.KeySerial1');
    elseif (strcmp(ctype,'glnx86') || strcmp(ctype,'glnxa64'))
        % For Linux machine with native serial port on motherboard. Just call the
        % serial command with an invalid param to get the list of possible serial
        % ports (/dev/ttyS*)
        display('Linux detected - opening serial port /dev/ttyS0');
        out = serial('/dev/ttyS0');
    else
        error('Computer architecture not recognized - default serial port options only specified on Mac and Linux');
    end
else
    % serial port string specified in function call
    out = portstring;
end

fopen(out);

w = Screen('OpenWindow',0);
LoadIdentityClut(w);

Screen(w,'FillRect',0);

Screen('Flip',w);

calibVals = cell(length(vals),1);

for i = 1:3;
    pause(1);
    beep
end

pad = zeros(1,length(vals));
if (strcmp(color,'gray'))
    rgbvals = [vals;vals;vals];
elseif (strcmp(color,'red'))
    rgbvals = [vals;pad;pad];
elseif (strcmp(color,'green'))
    rgbvals = [pad;vals;pad];
elseif (strcmp(color,'blue'))
    rgbvals = [pad;pad;vals];
elseif (strcmp(color,'other'))
    rgbvals = vals;
    if (size(rgbvals,1) ~= 3)
        error('vals must be a 3xN list of RGB values')
    end
end

for i = 1:length(vals)
    if (nargin < 3) % no size specified in function call
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
       
    pause(5);
    
    fprintf(out,'r');
    calibVals{i} = fgetl(out);
end

delete(instrfind);

sca

calibVals = cellfun(@pullOutNumber,calibVals);

end

function n = pullOutNumber(s)

s = s(3:strfind(s,'cd/m2')-1);

switch s(end)
    case ' '
        n = str2double(s);
    case 'm'
        n = str2double(s(1:end-1))/1000;
    case 'k'
        n = str2double(s(1:end-1))*1000;
    case 'u' % this one is likely wrong
        n = str2double(s(1:end-1))/10000;
    otherwise
        disp('Problem');
end

end
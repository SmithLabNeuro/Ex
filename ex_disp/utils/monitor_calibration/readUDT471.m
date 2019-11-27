%function lumVals = readUDT471(grayval, nvals, interval,portstring)
%
% This function shows the "grayval" and then for nvals readings with
% interval seconds between them it reads the value on the photometer.
% Useful for measuring the time course over which the CRT warms up.
%
% The serial port string for opening (e.g., '/dev/ttyS0'). If no string is
% specified, it assumes '/dev/cu.KeySerial1' on the mac and '/dev/ttyS0/'
% on Linux.
%
% Example Syntax:
% readUDT471(255,100,5);
%
% Matthew Smith
% Revised: 20120413

function lumVals = readUDT471(grayval, nvals, interval,portstring)

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

lumVals = cell(length(nvals),1);

for i = 1:3;
    pause(0.25);
    beep
end

for I=1:nvals
    Screen(w,'FillRect',grayval);
    Screen('DrawText',w,num2str(I),100,100);
    Screen('Flip',w);
       
    pause(interval);
    
    fprintf(out,'r');
    lumVals{I} = fgetl(out)
end

delete(instrfind);

sca

lumVals = cellfun(@pullOutNumber,lumVals);

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
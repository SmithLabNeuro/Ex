function microStim(t)
%
% Sends a TTL pulse out the "microStim" BNC line from the breakout box
% 
% Note: This is a backwards-compatible function from the November
% 2015 Linux upgrade for Ex. If you're on a PC it will merely use
% the old mex code (now renamed "winMicroStim"). If you're not on a
% PC, it will use the new comedi-based code.
%
% Modified 11Aug2017 by Matt Smith to remove windows, linux only now
%
% OLD WINDOWS COMMENTS HERE: 
% outputs a digital pulse on bit 3 of the AUXPORT. If 't' is specified, the
% duration of the pulse is 't' ms, otherwise it is 10ms.
%
% recompile with this command:
% mex winMicroStim.c cbw32.lib

global params;

% new comedi-based Linux digital output function
unixSendPulse(params.digOut1,t);

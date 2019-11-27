function [ output_args ] = timerDebugError( src,evt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global LE
LE{1} = lasterror;
LE{2} = evt;
LE{3} = src;

end


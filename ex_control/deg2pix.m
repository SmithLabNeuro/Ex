function pix = deg2pix(deg,scrd,pixpercm)
%DEG2PIX takes a number of degrees with the screen distance (in cm)
% and the pixels per cm and returns the pixels on the screen
%
% if only the degrees are passed in, it uses the globals to find the last two
% values.
%
% deg2pix(deg)
% deg2pix(deg,scrd,pixpercm)

% Matthew A. Smith
% Revised: 20130805

global params;

if isempty(params)
    exGlobals;
end
    
if (nargin == 1)
    scrd = params.screenDistance;
    pixpercm = params.pixPerCM;
end

angle = (pi/180) * deg;

pix = tan(angle) * scrd * pixpercm;

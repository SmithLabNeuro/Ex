function dva = pix2deg(pix,scrd,pixpercm)
%PIX2DEG takes a number of pixels with the screen distance (in cm)
% and the pixels per cm and returns the degrees of visual angle
%
% if only the pixels are passed in, it uses the globals to find the last two
% values.
%
% pix2deg(pix)
% pix2deg(pix,scrd,pixpercm)

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

d = pix./pixpercm;
angle = atan(d./scrd);

dva = (180/pi) * angle;
function [xVal, yVal, zValAng, buttonPress] = sampleHallEffectJoystick()

global params behav

buttonAnalogChannel = 6; % follow colors on joystick pinout
xAnalogChannel = 3; % follow colors on joystick pinout
yAnalogChannel = 4; % follow colors on joystick pinout
zAnalogChannel = 5; % follow colors on joystick pinout

button = unixGetAnalogInput(buttonAnalogChannel); % in mV
xVal = unixGetAnalogInput(xAnalogChannel); % in mV
yVal = unixGetAnalogInput(yAnalogChannel); % in mV
zVal = unixGetAnalogInput(zAnalogChannel); % in mV

if ~isfield(behav, 'zValAngle')
    behav.zValAngle = nan(100000, 1);
end

mvPerVolt = 1000;
posBaselineVolt = 2.5;
% xVal & yVal are converted into unit-less values that range from -1 to 1.
xVal = (xVal/mvPerVolt-posBaselineVolt)/posBaselineVolt;
yVal = (yVal/mvPerVolt-posBaselineVolt)/posBaselineVolt; 

mvPer45Degrees = 2500;
zValAng = (zVal-params.hallEffectZBaseline)/mvPer45Degrees*45; % convert to degrees turn

behav.zValAngle(find(isnan(behav.zValAngle), 1)) = zValAng;

% the button (at least on our original joystick) is very finicky, so I'm
% just defaulting to it not being pressed... can return later once we
% figure out any electrical problems (or if in a new joystick this
% works...)
buttonPress = false;

% if button > 0 && button < 20
%     buttonPress = true;
% else
%     buttonPress = false;
% end
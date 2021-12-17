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

if button < 10
    buttonPress = true;
else
    buttonPress = false;
end
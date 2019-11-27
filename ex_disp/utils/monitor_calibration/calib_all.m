% Calibrate gray, red, blue, and green. Make sure to warm up the monitor
% first (minimum 1 hr). Then output is saved to a .mat file as well as
% a .txt file for use with showex. For now, this script collects the
% individual color values, but then just uses gray for calibration.

clear all
calibvals = 0:255;
lv = calibrateMonitor(calibvals,'gray');
rv = calibrateMonitor(calibvals,'red');
bv = calibrateMonitor(calibvals,'blue');
gv = calibrateMonitor(calibvals,'green');

suffix = datestr(now,'yyyy.mmm.dd.HH.MM.SS');
calibdata = [calibvals' gv];

% save .mat file
save(['lumvals_',suffix,'.mat']);

% save photometervals file
save(['photometervals_',suffix,'.txt'],'calibdata','-ascii');


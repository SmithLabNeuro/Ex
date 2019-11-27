disp('The Change History for Ex and Ex-related functions are stored in comments in this file');
% ChangeLog.m
%
% 2016/05/16 by Matt Smith - changed unixGetEyes to return the voltage in
% mV, that way the value is the same as Ripple reads. Also changed
% wins.voltScale in runex to 5000, which was necessary given that the
% voltage range was different
%
% 2016/03/18 by Matt Smith - released the new version of RunEx/ShowEx. Many
% changes:
% (1) ShowEx now behaves differently with msgAndWait - it waits for the flip
% before returning the message. This means the codes will be better aligned
% with the diode.
% (2) Eye calibration data is now stored in raw volts. This makes
% recreating the calibration easier and more direct. "wins" is no longer
% sent via sendStruct, because it contains only variables relevant to the
% local control display. pixelsPerMV, for instance, which used to be needed
% to recreate eye position, is gone.
% (3) all "wait" functions are now optimized to check eye position the same
% way as rapidly as possible. 
% (4) 
%
% 2016/02/26 by Matt Smith - changed trialData display to show mouse mode
% key and resolution. Added sending of current condition and trial number
% to the ex_ function via 'e' struct.
%
% 2015/12/-- by Matt Smith - major update to core functions to support
% linux, add use of Priority call, change timer functions (remove
% plotEyes/plotMouse), reroute all position checks via the samp function,
% and more. Also all traces of the hist-related code (for user display of
% spike histograms) are removed, and other things are cleaned up (such as
% creating the local exData directory, etc) and made cross-platform
% compatible. Now uses the comedi library on linux and the cbw32 library on
% windows (and cbw64 for 64-bit Windows 7 compatibility). In the
% end, windows still needs work if we want to support it.
%
% 2015/12/04 by Adam Snyder - (1) removed nonlinear (polynomial and
% spherical) calibration capabilities.
%
% 2015/07/14 by Adam Snyder - (1) added reinitialization of xippmex if
% needed after matlapUDP call on close of runex. This should fix a problem
% where microstimulation fails when restarting runex after initializing.
%
% 2014/12/10 by Adam Snyder - (1) added a conditional to use
% experiment-specific 'juiceX' values so that the reward schedule can be
% easily randomized/parameterized.
%
% 2014/08/07 by Adam Snyder - (1) fixed bug where trial counter was
% incorrect when params.conditionFrequency was non-uniform.
%
% 2014/01/05 by Adam Snyder - (1) changed specification of calibration file
% to a full path.
%
% 2013/11/06 by Adam Snyder - (1) disabled automatic calibration.
%
% 2013/11/04 by Adam Snyder - (1) added support for the calibration dot
% size as a parameter (e.g., in the rig- or subject-specific xml option
% file). (2) added support for automatic calibration by pressing 'a' at the
% main screen.
%
% 2013/10/29 by Adam Snyder - (1) incorporated changes to calibration
% routine to enable polynomial regressors on x and y.
%
% 2013/09/16 by Adam Snyder - (1) changed calibration files so that they
% are stored locally on the control computer. The mouse mode calibration is
% loaded by default if there is not a calibration file saved locally yet.
%
% 2013/09/12 by Adam Snyder - (1) added support for rig-specific 'globals'.
% (2) fixed bugs I had introduced into the calibration routine. (3) changed
% character-matching conditionals to switch/case statements, arranged cases
% alphabetically. (4) moved the experimental control stuff to a
% subfunction, which makes the code a little easier to navigate. (5)
% stopped timer functions before changing them during mouse toggle, then
% restart them. This seems to have helped the intermittant freezing when
% toggling. (6) Made automatic local directory handling, and added a
% command that sets local directory stuff on the showex side.
%
% 2013/09/12 by Adam Snyder - (1) moved wins settings from exGlobals.m to
% runex.m, removed references to histogram window. Now detects the Ex
% control monitor resolution and resizes the display layout accordingly.
% (2) hitting 'q' within a block now "pauses"; pressing 's' again continues
% the block from the next trial (the old behavior had been to start the
% current block over again).
%
% 2013/09/09 by Adam Snyder - (1) changed error handling when running the
% experiment. Rather than just having a try/catch around the ex-file
% execution, the try is now around everything that happens when getChar is
% 's'. Error handling was moved to the new subfunction runexError, which
% prints the error message on the Ex Control screen, and returns keyboard
% control to the user to exit gracefully. The full error message and stack
% is printed in the command window. (2) Any numeric digit keypress now
% adjusts the number of juice clicks (juiceX). It had been just the numbers
% 1-7. (3) removed miscellaneous calls to histogram-related functions.
%
% 2013/09/05 by Adam Snyder - (1) added support for subject XML files.
% Subject-specific settings will overwrite the default settings in the
% template XML files. You only need to specify subject settings that you
% want to differ from the defaults. You can also specify subject-specific
% globals. This change involved making massive changes to readExperiment,
% and the creation of two new functions: readSubject and parseExperiment.
% The framework exists to add arbitrary layers of priority for settings
% (e.g., subject>default or subject>rig-specific>default, etc.), although
% some little changes will be needed to implement that (which I plan to do
% shortly). (2) removed directory checking. Directories should be on the
% MATLAB search path now. May need to check where datafiles are being
% saved, though.
%
% 2013/09/03 by Adam Snyder - (1) calibration can now be set to present the
% calibration points in a random order by setting
% params.randomizeCalibration in the 'exGlobals' file to logical true.
% Calibration is not randomized by default. (2) moved toggling between
% mouse and monkey mode to a subfunction --will come back to fix this soon.
% (3) changed 'drawCalibration' so that if no input is passed it draws dots
% for all the points. (4) added support for trial aborts when showex misses
% a message, which would otherwise cause message and wait to hang. This
% involved changing msgAndWait, too.
%
% 2013/08/13 by Adam Snyder - (1) prompt for new subject ID if current value in
% exGlobals.m isn't correct (you can still ctrl+c to break and change globals if
% desired). (2) Miscellaneous code cleanup. (3) Rethrow the specific error
% when reading an XML file when that problem is encountered. (4) added code
% cells to improve navigation. Reorganized code to fit the cell
% categorizations. (5) close open udp sessions before trying to initialize
% a new udp session. (6) now save the name of the control machine with the
% data (as a parameter). (7) moved calibration to a subfunction.
%
% 2012/10/22 by Matt Smith - enables support for automatic output
% file naming. Also fixes the path issues.
%
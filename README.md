—— <BR>
<B>Ex - Smith Lab Experimental Software</B><BR>

This repository contains software for online experimental control of neurophysiology experiments. It is primarily maintained and used by the Smith Laboratory (<A HREF="www.smithlab.net">www.smithlab.net</A>), but it is also in use by several other labs. Ex leverages the Psychophysics Toolbox (<A HREF="psychtoolbox.org">psychtoolbox.org</A>) for generating visual displays, and is written in Matlab (<A HREF="www.mathworks.com">www.mathworks.com</A>) native code with some C-code compiled in Matlab for speed.<P>

—— <BR>
<B>Ex History</B> <BR>

Ex was created by Ryan Kelly and Matt Smith in ~2010, with some building blocks started a little before that. Ryan began this work as a PhD student in Tai Sing Lee’s lab and wrote most of the initial version of Ex, and then when he graduated Matt took over sole responsibility for code development and maintenance.<P>

Many people have contributed to Ex over the years, but a few notables include:<BR>
Adam Snyder - extensive improvements throughout the system as it matured from its initial version<BR>
Ryan Williamson - implemented the initial version of BCI code<BR>
Xiaohan Zhong - improvements throughout, in particular in XML file reading and system initialization<BR>
Emilio Salazar-Gatzimas - formalized and extended the BCI realtime framework, implemented an experimental notes database using SQLite, as well as other improvements<BR>

——<BR>
<B>Ex Overview</B><BR>

Ex involves, at a minimum, two computers, a “control” computer and a “graphics/display” computer (each running Linux). These two computers communicate via a direct ethernet cable (using UDP). The control computer is the one where the operator interacts the most, setting in motion the experiments and watching the eye movements and behavior of the subject. The graphics computer mostly just runs a single program that sits and waits for commands about what to show. There are optionally two other computers that can be in place, a data recording computer and a bci/real-time analysis computer. For these features, the Ex system is tightly integrated with software and hardware provided by Ripple (<A HREF="www.rippleneuro.com">www.rippleneuro.com</A>) for data acquisition.<P>

——<BR>
<B>Ex Repository Directory Structure</B><BR>

ex - Contains individual experimental control functions, named “ex_taskname.m”. Two examples are provided, ex_demotest.m and ex_timingtest.m<BR>
ex_bci - Code related to the BCI/realtime analysis, including the “rtex.m” main function<BR>
ex_control - Code related to the main experimental control computer, including the “runex.m” function<BR>
ex_data - Code for online data recording control, including the “recordex.m” function<BR>
ex_disp - Code for the display computer, including the “showex.m” function<BR>
stim - Library of functions for showing stimuli, for example “stim_grating.m”<BR>
xippmex - Code library provided by Ripple for interfacing with their hardware in Matlab<BR>
xml - XML files that provide parameters to run experiments (these pair with the “ex_taskname.m” functions)<BR>


NOTE: You can use the “ex” and “xml” directories to store your experimental control and configuration files, but best practice is to create an additional directory structure in a directory called “Ex_local” where you put those files. The location for that directory, and for some other key directories, is set in the “ex_control/exGlobals.m” file<P>

——<BR>
<B>Hardware requirements</B><BR>

At a minimum, two computers are required. We run Ubuntu Linux, usually the most recent stable release (labeled LTS), along with the most recent PsychToolbox and Matlab. Some attention has to be paid to make sure the OS and software are fully compatible with each other - sometimes new versions break this temporarily.<P>

For the control and display computers, a typical purchase would be a relatively powerful recent computer (Intel i7 chip) with at least 16GB of RAM and a good SSD, and two additional hardware specs: (1) A multiport ethernet card (a single additional port is sufficient for the display, but usually a 2-port or 4-port card is better for the control depending on how you’re going to configure your system); (2) A good graphics card (AMD Radeon cards in the Polaris family are recommended - e.g., AMD Radeon WX3200). In addition, the control computer needs an analog & digital I/O card (Measurement Computing PCIe-DAS1602/16 board). You’ll also need an ethernet cable to connect the extra port from these two computers directly.<P>

——<BR>
<B>Getting started</B><BR>

The first step is configuring the two computers. In brief, this process typically involves:<BR>
(1) Install Ubuntu Linux and configure with a low latency kernel<BR>
(2) Install all other key software packages (Matlab, Ex, Psychtoolbox via NeuroDebian, comedi for I/O card control, samba for file sharing, git for repository access)<BR>
(3) Configure Ubuntu (disable screensaver and automatic suspend, turn off notifications, turn off automatic updates, set refresh rate and resolution to desired values)<BR>
(4) Configure Matlab (add the appropriate Ex directory [ex_control, ex_disp, etc] to path, run “PsychLinuxConfiguration” from PsychToolbox)<BR>
(5) Configure ethernet ports on each computer for direct communication (use 192.168.1.10 for control and 192.168.1.11 for display by default)<BR>
(6) copy the exlocal_template files to an Ex_local directory for your computer ('cp -r misc/exlocal_template ../Ex_local')

Once the ethernet ports of the two computers are configured properly (use “ping” at the terminal to test the connection), you should first run “showex” on the display computer, followed by “runex” on the control computer. The control computer checks for communication with the display before proceeding. Example syntax for the “runex” command would be as follows:<P>

runex(‘demotest.xml’);<BR>

When you first execute this “runex” command an “exData” directory will be created in the user’s root directory, where Matlab-native data files (.mat) with a record of the task parameters and behavioral events will be saved (one file for each time you execute “runex”). There are other optional arguments to runex, see the comments at the top of runex.m to learn more about them.<P>

——<BR>
<B>Controlling runex</B><BR>

When runex successfully loads and communicates with showex, it will bring up an interface. Keyboard commands indicate your control options. A simple first step is to press ‘m’ to enter “mouse mode”, where the eye position for the task is controlled by the operator’s mouse. Then press ’s’ to start the task. The user interface indicates the valid keystroke options to control runex behavior during the task, and to halt/change/restart the task, and the moving cursors at the top show the raw eye signal input (on the left) and the calibrated eye signal (on the right).<P>

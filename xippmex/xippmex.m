% xippmex_help.m Help file for XIPPMEX MEX-file
%
% XIPPMEX Matlab interface to NIP and Trellis software through XIPP.
%
% usage: xippmex(cmdstr [, args])
%                 'time' - display latest NIP time
%                 'spike' - get recent spike counts and times
%                 'cont' - get continuous data
%                 'close' - close UDP socket and delete cached data
%                 'stim' - send stim control string
%                 'signal' - enable or disable signals
%                 'digin' - retrieve digital inputs
%                 'digout' - control digital outputs
%                 'filter' - modify and retrieve filter information
%                 'stimseq' - complex control of stimulation
%                 'priority' - change xippmex thread priority
%                 'opers' - find Trellis operators on the network
%                 'fastsettle' - control or display NIP fast settle
%                 'trial' - control file save on Trellis operators
%                 For more info, type command strings without arguments
%
% xippmex version 1.2.3.24
%
%  MEX-File function

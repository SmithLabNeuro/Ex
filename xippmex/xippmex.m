% xippmex_help.m Help file for XIPPMEX MEX-file
%
% XIPPMEX Matlab interface to NIP and Trellis software through XIPP.
%
% usage: xippmex(cmdstr [, args])
%                 'time' - display latest NIP time
%                 'elec' - retrieve electrode numbers by type of headstage
%                 'spike' - get recent spike counts and times
%                 'spike-thresh' - set and receive spike thresholds
%                 'cont' - get continuous data
%                 'digin' - retrieve digital inputs
%                 'digout' - control digital outputs
%                 'stim' - send stim control string
%                 'stimseq' - complex control of stimulation
%                 'signal' - enable or disable signals
%                 'filter' - modify and retrieve filter information
%                 'fastsettle' - control or display NIP fast settle
%                 'lowcorner' - set an electrode's hardware filter low corner
%                 'adc2phys' - set an electrode's ADC resolution

%                 'trial' - control file save on Trellis operators
%                 'impedance' - trigger impedance measurement
%                 'close' - close UDP socket and delete cached data
%                 For more info, type command strings without arguments
%
% xippmex version 1.5.0.91
%
%  MEX-File function

function ordering = createOrdering(expt,varargin)
%createOrdering: a function to create ordering vector for Ex trial
%presentation. -Adam C. Snyder (adam@adamcsnyder.com)
global behav
%% Parse inputs
p = inputParser;
p.addRequired('expt',@iscell);                    %'expt' is a cell array returned by readExperiment
p.addParamValue('blockRandomize',true);           %whether or not to randomly permute the conditions
p.addParamValue('conditionFrequency','uniform');  %frequency of conditions, use a 1xNumConditions vector of integers to specify frequency of each condition within a block
p.addParamValue('numBlocksPerRandomization',1);   %randomize over 1 or more blocks (i.e., if numBlocksPerRandomization ==2, then concatenate two blocks before shuffling, etc.)
p.addParamValue('exFileControl','no');            %whether the ordering list is completely under EX file control or not. Specify 'no' or 'yes', or an integer which specifies which condition to present for the first trial (a random first trial is selected for 'yes')
p.parse(expt,varargin{:});
%% EX file control
if ~strcmpi(p.Results.exFileControl,'no')
    if ~isfield(behav,'ordering') %if ordering is not specified in the behav struct by the EX file (e.g., for the very first trial)
        if isnumeric(p.Results.exFileControl)
            behav.ordering = p.Results.exFileControl; %use the starting condition(s) specified           
        else
            behav.ordering = randi(numel(expt)); %If no starting condition is specified, then just pick one at random
        end;
    end;
    ordering = behav.ordering; %just use behav.ordering if it is provided
    return; %done: exFileControl trumps everything else
end;
%% Block randomization
if isnumeric(p.Results.conditionFrequency)
    assert(numel(p.Results.conditionFrequency)==numel(expt),'RUNEX:ConditionFrequencyMismatch','If specifying condition frequencies, the number of frequencies specified must match the number of conditions');
    ordering = zeros(1,sum(p.Results.conditionFrequency)); %preallocate
    frequencies = [0 p.Results.conditionFrequency]; %preappend a zero to make the cumsum thing work without an exception for the first condition...
    for fx = 2:numel(frequencies)
        startInd = sum(frequencies(1:fx-1))+1; endInd = startInd+frequencies(fx)-1;
        ordering(startInd:endInd) = fx-1;
    end;
else
    ordering = 1:numel(expt); %just use uniform single-instance-per-block frequency by default
end;
ordering = repmat(ordering,1,p.Results.numBlocksPerRandomization); %make the requested number of copies of blocks over which to randomize
if p.Results.blockRandomize, ordering = ordering(randperm(numel(ordering))); end; %if block randomization is false, then just return the conditions in order I guess...
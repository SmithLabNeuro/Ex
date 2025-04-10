function decoderFileLocationAndName = ...
    calibrateKalmanVelocityDemoDecoderForBci...
    (~, nevFilebase, nevFilesForTrain, trainParams, subject)

% This script calibrates a velocity Kalman filter (vKF). It implements the
% followinf steps:
% Step 1:  Read and organize the data
% Step 2:  Preprocess neural data. choose releevnt traisl, aligning the data
%          to BCI onset, binning
% Step 3:  Choose channels with sufficient quality (based on firing rate, 
%          fano factor, and coincident spikes)
% Step 4:  Calculate spike count matrix 
% Step 5:  Calculate assumed velocity for vKF fitting
% Step 6:  Z-score neurons and fit FA model to dimensionality reduce neural 
%          activity
% Step 7:  Fit Kalfman filter
% Step 8:  Converge Kalman Gain
% Step 9:  Calculate M matrices
% Step 10: Calculate trajectories based on current model & plot
% Step 11: Save model

global params codes

VEL_DIM = 2; % dimensionality of the cursor velocity (horizontal and vertical)
alignSpikesToCode = trainParams.bciStartCode; % code to start BCI control
timeBefore = trainParams.timeBeforeMs / 1000; %sec, time to take before code
timeAfter = trainParams.timeAfterMs / 1000; %sec, time to take after code
binSize = trainParams.binSizeMs / 1000; %sec
numBins = floor((timeAfter-timeBefore)/binSize); % number of bins these paramenetes imply

%% Step 1 - READ DATA

% Load Nev Files
firstTimeptShift = 0;
spikesInit = zeros(0, 3);
wavesInit = zeros(52, 0);
for i = 1:length(nevFilesForTrain)
    
    currNevFile = nevFilesForTrain{i};
    fprintf('Reading NEV file:  %s \n', currNevFile)
    [spikes, waves] = readNEV(currNevFile);
    
    % apply timeshift for multiple NEV files
    spikes(:, 3) = spikes(:, 3) + firstTimeptShift;
    firstTimeptShift = spikes(end, 3) + 1; % shift of 1s
    spikesInit = cat(1, spikesInit, spikes);
    wavesInit = cat(2, wavesInit, waves);

end

fullNevDat = {spikesInit, wavesInit};

[~,nevLabelledData] = runNASNet(fullNevDat,...
    trainParams.gamma,...
    'channels', trainParams.rippleChannelNumbersInBci,...
    'netFolder', params.nasNetFolderDataComputer,...
    'netname', trainParams.nasNetwork);

dat = nev2dat(nevLabelledData, 'nevreadflag', true);

dat = bcidemo_getData(dat);

%% Step 2: Preprocess neural data

%% A - by rate anf fano factor
relevantResultCodes = cellfun(@(x) codes.(x), trainParams.trainingResultCodes);
inxReleventTrials = find(ismember([dat.trials.outcome],relevantResultCodes));
dat.trials = dat.trials(inxReleventTrials);

rates = nan(1,length(trainParams.rippleChannelNumbersInBci));
fanoFactor = nan(1,length(trainParams.rippleChannelNumbersInBci));

for i = 1:length(trainParams.rippleChannelNumbersInBci)

    rates(i) = bcidemo_getChanRate(dat,trainParams.rippleChannelNumbersInBci(i));

    fanoFactor(i) = bcidemo_getChanFF(dat,trainParams.rippleChannelNumbersInBci(i));
    
end

disp(['Number removed due to rate = ' num2str(sum(~(rates>trainParams.firingRateThreshold)))])

disp(['Number removed due to FF = ' num2str(sum(~(fanoFactor < trainParams.fanoFactorThreshold)))])

channelsKeep = trainParams.rippleChannelNumbersInBci(find(...
    rates>trainParams.firingRateThreshold & ...
    fanoFactor < trainParams.fanoFactorThreshold));


%% B - remaining channels - chekced for coinsident spikes

disp('Calcuting coincidence matrix')
coincidenceSpiking = nan(length(channelsKeep),length(channelsKeep));

% Note: this will not return a symatrical matrix because while the number
% of coincident spikes between to channels is equal the total number of
% spikes per channel will be different and so the returned fraction will be
% different. 
% The (i,j)th element of the matrix is the fraction of spikes in channel i
% that have a concurrent spike in channel j. 

for i = 1:length(channelsKeep)

    coincidenceSpiking(i,:) = bcidemo_getCCSpiking (dat, channelsKeep(i),...
        channelsKeep, trainParams.coincTimeMs/1000);
    coincidenceSpiking(i,i) = nan;
end

disp('Removing coincidenct channels')

remove_channels = [];

while any(coincidenceSpiking>trainParams.coincThresh,"all")

    [max_col,i_mx_col] = max(sum((coincidenceSpiking > trainParams.coincThresh),2));
    [max_row,i_mx_row] = max(sum((coincidenceSpiking > trainParams.coincThresh),1));
    if max_col>max_row
        i_mx = i_mx_col;
    else 
        i_mx = i_mx_row;
    end

    remove_channels = [remove_channels,channelsKeep(i_mx)];
    coincidenceSpiking(:,i_mx) = nan; 
    coincidenceSpiking(i_mx,:) = nan; 
end

disp(['Number removed due to CC = ' num2str(length(remove_channels))])

channelsKeep = setdiff(channelsKeep,remove_channels);

disp("Kept " + num2str(length(channelsKeep) + " channels for decoder"))

%% Step 4: Calculate spike count matrix

fr = nan(length(channelsKeep),length(dat.trials),numBins);

for i = 1:length(channelsKeep)

    raster = bcidemo_calculateRasters(dat, alignSpikesToCode, timeBefore, timeAfter, ...
        1:length(dat.trials), channelsKeep(i));
    binned = bcidemo_binRaster(raster,binSize,...
        timeBefore,...
        timeAfter);
    fr(i,:,:) = binned';
end

spikeCounts = fr * binSize;



%% Step 5: Calculate assumed velocity

% assumptions about speed
assumedCursorSpeed = trainParams.assumedCursorSpeed; % px per sec 
assumedVel = nan(VEL_DIM, length(dat.trials),numBins);

for i = 1:length(dat.trials)

    targetAngle = dat.trials(i).target;

    assumedVel(1,i,:) = cosd(targetAngle) * ones(1,numBins) * assumedCursorSpeed;
    assumedVel(2,i,:) = sind(targetAngle) * ones(1,numBins) * assumedCursorSpeed;

end

% figure; hold on
% for i = 1:length(dat.trials)
%     plot(squeeze(cumsum(assumedVel(1,i,:))) * binSize,...
%         squeeze(cumsum(assumedVel(2,i,:))) * binSize,...
%         '-o', 'DisplayName', num2str(dat.trials(i).target))
% end
% legend show

assumedVel = reshape(assumedVel,[VEL_DIM,length(dat.trials)*numBins]);

%% Step 6: Z-score and fit FA model

X = reshape(spikeCounts,[length(channelsKeep),length(dat.trials)*numBins]);

X_mean = mean(X,2);
X_std = std(X,0,2);
X = zscore(X,0,2);

[estFAParams, ~] = fastfa(X,trainParams.numberFaLatents);
[lat, ~, beta] = fastfa_estep(X, estFAParams);
lat = lat.mean;

%% Step 7: Kalman Model Fit

T = size(lat,2);

A = eye(VEL_DIM);

Q = eye(VEL_DIM) * trainParams.kalmanQ;

C =  (lat*assumedVel')* inv(assumedVel*assumedVel');

d = (1/T) * (sum(lat,2) - C * sum(assumedVel,2));

R = (1/T) * (lat-d-C*assumedVel) * (lat-d-C*assumedVel)';

%% Step 8: Kalman Gain

NUM_ITERATION_STEPS = 100;
TOL = 10^-10;
Sigma_t_t = cov(assumedVel');
convergenceCheck = nan(NUM_ITERATION_STEPS,VEL_DIM*trainParams.numberFaLatents);
kalmanGain = randn(size(2,10));

for i = 1:NUM_ITERATION_STEPS

    Sigma_t_tprev = A*Sigma_t_t*A'+ Q;
    
    convergenceCheck(i,:) = kalmanGain(:);

    kalmanGain = Sigma_t_tprev * C' * ...
    inv(C * Sigma_t_tprev * C' + R);
    

    Sigma_t_t = Sigma_t_tprev - kalmanGain * C * Sigma_t_tprev;

end

if any(abs(convergenceCheck(end,:) - convergenceCheck(end-1,:)) > TOL)
    disp('Kalman Gain not converged !!')
else
    disp('Kalman Gain converged')
end
 
%% Step 9: M matrices

M0 = - kalmanGain * d;
M1 = A - kalmanGain * C * A;
M2 = kalmanGain*beta;


%% Calculate trajectories based on current model & plot

figure; hold on

directions = [dat.trials.target];
unique_directions = unique(directions);
INITAL_VEL = [0,0]';
added_to_legend = false(1, length(unique_directions)); % Keep track of legend additions
colors = varycolor(length(unique_directions));

for j = 1 : length(unique_directions)

    inx = find(directions == unique_directions(j));
    leg{j} = num2str( unique_directions(j));
    
    for i = inx

        trajectory = zeros(VEL_DIM,1+size(spikeCounts,3));
        currVel = INITAL_VEL;

        for t = 1:size(spikeCounts,3)

            currVel = M0 + M1*currVel+ M2*((squeeze(spikeCounts(:,i,t))-X_mean)./X_std);
            trajectory(:,t+1) = currVel;

        end

        x = [cumsum(trajectory(1, 1:end))]' * binSize;
        y = [cumsum(trajectory(2, 1:end))]' * binSize;
        % Only add the first instance of the direction to the legend
        if ~added_to_legend(j)
            plot(x, y, 'Color', colors(j, :), 'Marker', 'o', 'MarkerSize', 10, 'DisplayName', leg{j});
            added_to_legend(j) = true;
        else
            plot(x, y, 'Color', colors(j, :), 'Marker', 'o', 'MarkerSize', 10, 'HandleVisibility', 'off');
        end
        
        
    end
end

axis([-300 300 -300 300])
legend(leg)
%% save model
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(params.bciDecoderBasePathDataComputer, bciDecoderRelativeSaveFolder);

subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
bciDecoderSaveName = sprintf('%s%sKalmanBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));

zScoreSpikesMat = true;
zScoreLatentMat = false;
estParams = estFAParams;
K = kalmanGain;
gamma = trainParams.gamma;
nasNetName = trainParams.nasNetwork;

save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName),...
    'M0', 'M1', 'M2', ...
    'channelsKeep',...
    'A', 'Q', 'C', 'R', 'beta',...
    'K',...
    'zScoreSpikesMat', 'zScoreLatentMat',...
    'estParams', 'nevFilebase',...
    'nevFilesForTrain',...
    'nasNetName', 'gamma',...
    'X_mean', 'X_std',...
    'trainParams');

decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);



 
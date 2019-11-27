function table = makeGammaTable(rgbvals,luminance,nvals,plotflag)
% function table = makeGammaTable(rgbvals,luminance,nvals,plotflag)
% 
% Makes a gamma correction table from luminance readings using a
% power law function fit. Can be modified to use the FitGamma
% function in PsychToolbox if desired.
%
% rgbvals is a set of input values for the color of the display
% luminance is the readings from a photometer for the rgbvals
% nvals is the number of entries in the LUT
% plotflag 1 for plots, 0 for no plots
%
% The output, table, is a 3-column (r,g,b) set of values ranging
% from 0-1 and has a default of 256 rows.
%  
% In psychtoolbox, the table is loaded this way:
%  Screen('LoadNormalizedGammaTable', ScreenNumber, table, [loadOnNextFlip], [physicalDisplay]);
%
%   Matthew A. Smith
%   Revision: 20080828
%
% 2016mar16 - suppress text output from lsqcurvefit - MAS
%

if (min(size(rgbvals)) > 1)
    warning('Monochrome calibration only!');
    return
end

if (nargin < 4)
    plotflag = 0;
end

if (nargin < 3)
    nvals = max(max(rgbvals))+1;
end

maxIndex = nvals-1; % +1 = total num of indices
indexNormal = [0:.001:1]'; 
maxLum = max(max(luminance));

% Fit power function to the luminance data
xb = [0.1]; % starting value
lowerbound = 1e-5;
upperbound = Inf;
opts = optimset('Display','off');
[fitval resval] = lsqcurvefit(@powerfunc_makeGammaTable,xb,rgbvals/maxIndex,luminance/maxLum,lowerbound,upperbound,opts);
% interpolate the values of the fitted curve
out1 = powerfunc_makeGammaTable(fitval,indexNormal);

% Alternatively, we could use the psychtoolbox fitting function:
%
%if (exist('FitGamma') == 2)
%    % You can select a fitting method with the last argument 
%    % see help for FitGamma in Psychophysics Toolbox 
%    [out1 fitval message1] = FitGamma(rgbvals/maxIndex,  luminance/maxLum, indexNormal, 1);
%else
%    warning('PsychToolbox function FitGamma.m not available, defaulting');
%    return
%end

% Create the gamma correction table [compute the inverse function numerically]
orgTable = [maxIndex*indexNormal  maxLum*out1];
correctedTable = [];
for i = 0:maxIndex
    eqLum = i * maxLum/maxIndex; 
    numTable = max(find(orgTable(:,2) <= eqLum)); 
    correctedTable = [correctedTable ;  i round(orgTable(numTable,1))]; 
end

% now make the table with a range from 0-1
col = correctedTable(:,2)./max(correctedTable(:,2));
table = [col col col];

% plots for confirmation of table linearity
if (plotflag)
    % original luminance readings
    figure(1); 
    subplot(1,3,1);
    plot(maxIndex*[0:0.001:1]', maxLum*out1); box off;
    hold on; 
    scatter(rgbvals,luminance);
    ylabel('luminance');
    
    % confirm linearity 
    subplot(1,3,2);
    linearTable = []; 
    for i=0:maxIndex
        numTable = max(find(orgTable(:,1) <= correctedTable(i+1, 2))); 
        linearTable = [linearTable;  i  orgTable(numTable, 2) ]; 
    end 
    plot(linearTable(:,1), linearTable(:,2)'); box off;
    xlabel('table entry');
    
    subplot(1,3,3);
    plot(table(:,1)); box off;
end

%------------------------------------------------------------------------------
%%% BEGINNING OF INSERTED FUNCTION "powerfunc_makeGammaTable"

function [y]=powerfunc_makeGammaTable(x,xdata)
powerVal = x(1);
xdataMax = max(xdata);
y = (xdata/xdataMax).^powerVal;

%%% END OF INSERTED FUNCTION "powerfunc_makeGammaTable"
%------------------------------------------------------------------------------

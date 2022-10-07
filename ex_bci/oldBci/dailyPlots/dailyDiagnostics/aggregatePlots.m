function [  ] = aggregatePlots( startDate, endDate, savePath, varargin )
%AGGREGATEPLOTS Summary of this function goes here
%   Detailed explanation goes here

    subjectName = 'Pe';
    filePath = '/Volumes/Queen/data/pepelepew';
    savePlots = false;
    assignopts(who,varargin);

    % get all bcistruct file names
    fNames = dir(sprintf('%s/%s*bcistruct*.mat',filePath,subjectName));
    
    % figure out start and end dates
    startDT = datetime(startDate,'InputFormat','yyMMdd');
    endDT = datetime(endDate,'InputFormat','yyMMdd');
    allDT = startDT:endDT;
    
    % go through and get smoothed perf and acqTime
    smoothedAcq = struct('t',[],'vals',[]);
    smoothedCorr = struct('t',[],'vals',[]);
    currIdx = 1;
    for ii = 1:length(allDT)
        fprintf('Loading %s...\n',datestr(allDT(ii)));
        foundFile = false;
        for i_file = 1:length(fNames)
            currDT = datetime(fNames(i_file).name(3:8),'InputFormat','yyMMdd');
            if isequal(currDT,allDT(ii))
                foundFile = true;
                break;
            end
        end
        if ~foundFile
            fprintf('   Could not find file for this date. SKIPPING...\n');
            continue;
        end
        load(sprintf('%s/%s',filePath,fNames(i_file).name));
        fprintf('   Getting smoothed bci performance...\n');
        smoothedCorr(currIdx) = plotSmoothBCICorrect(bcistruct.bcidata);
        close(gcf);
        smoothedAcq(currIdx) = plotSmoothBCIAcquisitionTime(bcistruct.bcidata);
        close(gcf);
        currIdx = currIdx+1;
    end
    for ii = 1:length(smoothedCorr)
        smoothedCorr(ii).T = length(smoothedCorr(ii).t);
    end
        
    % compile results into matrices
    [~,maxIdx] = max([smoothedCorr.T]);
    maxT = 90;
    [perfMatrix,acqMatrix] = deal(nan(length(smoothedCorr),maxT));
    for ii = 1:length(smoothedCorr)
        currT = min(smoothedCorr(ii).T,maxT);
        perfMatrix(ii,1:currT) = smoothedCorr(ii).vals(1:currT);
        acqMatrix(ii,1:currT) = smoothedAcq(ii).vals(1:currT);
    end
    
    % plot results
    t = smoothedCorr(maxIdx).t;
    xIdx = floor(linspace(1,maxT,5));
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[0 1 1.5 2]);
    subplot(2,1,1);
    imagesc(perfMatrix,[0 .8]); colorbar
    set(gca,'XTick',xIdx,'XTickLabel',floor(t(xIdx)),...
        'YTick',[1 length(smoothedCorr)],'YTickLabel',datestr(allDT([1 length(allDT)])));
    xlabel('time (min)'); ylabel('session'); title('Percent correct');
    
    subplot(2,1,2);
    imagesc(acqMatrix,[0 1500]); colorbar
    set(gca,'XTick',xIdx,'XTickLabel',floor(t(xIdx)),...
        'YTick',[1 length(smoothedAcq)],'YTickLabel',datestr(allDT([1 length(allDT)])));
    xlabel('time (min)'); ylabel('session'); title('Acquisition time (ms)')
    
    if savePlots
        print([savePath '/aggregatePerf.png'],'-dpng');
    end
    
end


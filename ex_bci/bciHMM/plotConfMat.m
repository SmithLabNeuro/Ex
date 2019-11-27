function [  ] = plotConfMat( predStruct, whichType )
%PLOTCONFMAT Summary of this function goes here
%   Detailed explanation goes here

    y_true = predStruct.pred_both.y_true;
    uniqueY = unique(y_true);
    if strcmpi(whichType,'training')
        y_both = predStruct.pred_both.y_training;
        y_left = predStruct.pred_left.y_training;
        y_right = predStruct.pred_right.y_training;
    elseif strcmpi(whichType,'crossval')
        y_both = predStruct.pred_both.y_cv;
        y_left = predStruct.pred_left.y_cv;
        y_right = predStruct.pred_right.y_cv;
    else
        error('Invalid value for whichType ("training" or "crossval").')
    end
    
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 2]);
    
    % left
    subplot(2,2,1);
    [left_conf,angleOrder] = confusionmat(y_true,y_left);
    uniqueAngles = angleOrder(1:4:length(angleOrder));
    for ii = 1:length(angleOrder)
        left_conf(ii,:) = left_conf(ii,:) ./ sum(y_true==angleOrder(ii));
    end
    imagesc(left_conf,[0 1]); colorbar;
    xlabel('Predicted angle'); ylabel('True angle');
    set(gca,'XTick',find(ismember(angleOrder,uniqueAngles)),'YTick',find(ismember(angleOrder,uniqueAngles)));
    set(gca,'XTickLabel',uniqueAngles,'YTickLabel',uniqueAngles);
    title(sprintf('Left array, %s',whichType));
    
    % right
    subplot(2,2,2);
    [right_conf,angleOrder] = confusionmat(y_true,y_right);
    uniqueAngles = angleOrder(1:4:length(angleOrder));
    for ii = 1:length(angleOrder)
        right_conf(ii,:) = right_conf(ii,:) ./ sum(y_true==angleOrder(ii));
    end
    imagesc(right_conf,[0 1]); colorbar;
    xlabel('Predicted angle'); ylabel('True angle');
    set(gca,'XTick',find(ismember(angleOrder,uniqueAngles)),'YTick',find(ismember(angleOrder,uniqueAngles)));
    set(gca,'XTickLabel',uniqueAngles,'YTickLabel',uniqueAngles);
    title(sprintf('Right array, %s',whichType));
    
    % both
    subplot(2,2,3);
    [both_conf,angleOrder] = confusionmat(y_true,y_both);
    uniqueAngles = angleOrder(1:4:length(angleOrder));
    for ii = 1:length(angleOrder)
        both_conf(ii,:) = both_conf(ii,:) ./ sum(y_true==angleOrder(ii));
    end
    imagesc(both_conf,[0 1]); colorbar;
    xlabel('Predicted angle'); ylabel('True angle');
    set(gca,'XTick',find(ismember(angleOrder,uniqueAngles)),'YTick',find(ismember(angleOrder,uniqueAngles)));
    set(gca,'XTickLabel',uniqueAngles,'YTickLabel',uniqueAngles);
    title(sprintf('Both arrays, %s',whichType));
    
    % tuning curves
    subplot(2,2,4); hold on;
    extAngles = [uniqueY(end) uniqueY uniqueY(1)];
    corr_both = nan(length(y_true),1); corr_left = nan(length(y_true),1); corr_right = nan(length(y_true),1);
    for i_trial = 1:length(y_true)
        angIdx = find(uniqueY==y_true(i_trial));
        corr_vec = extAngles(angIdx:angIdx+2);
        corr_both(i_trial) = sum( ismember(corr_vec,y_both(i_trial)) ) > 0;
        corr_left(i_trial) = sum( ismember(corr_vec,y_left(i_trial)) ) > 0;
        corr_right(i_trial) = sum( ismember(corr_vec,y_right(i_trial)) ) > 0;
    end
    both_tuning = nan(length(uniqueY),1); left_tuning = nan(length(uniqueY),1); right_tuning = nan(length(uniqueY),1);
    for i_ang = 1:length(uniqueY)
        tmpIdx = y_true==uniqueY(i_ang);
        tmpCorrLeft = corr_left(tmpIdx); tmpCorrRight = corr_right(tmpIdx); tmpCorrBoth = corr_both(tmpIdx);
        left_tuning(i_ang) = sum(tmpCorrLeft)./length(tmpCorrLeft);
        right_tuning(i_ang) = sum(tmpCorrRight)./length(tmpCorrRight);
        both_tuning(i_ang) = sum(tmpCorrBoth)./length(tmpCorrBoth);
    end
    plot(uniqueY,left_tuning,'ro-');
    plot(uniqueY,right_tuning,'bo-');
    plot(uniqueY,both_tuning,'ko-','LineWidth',3);
    set(gca,'XTick',uniqueAngles); ylim([0 1]);
    xlabel('Target angle'); ylabel('Prediction accuracy');
    title(sprintf('%s pred acc (neighbor angles as correct)',whichType));
    legend('left','right','both','Location','Best'); legend boxoff;
    
end


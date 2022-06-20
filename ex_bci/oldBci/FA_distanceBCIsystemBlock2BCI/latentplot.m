function h = latentplot(dat,indices,latenttype,latentinds,plottype,maxtraj)

%% plotting code
%% Input variables
% indices must either be an array of doubles, a cell array of double arrays
% or a character array
% if indices is a character array
% 'b' = bcimissed
% 'g' = bcihit (bci hit)
% 'm' = missed
% 'h' = hit
% 'c' = calibration
% example with doubles indices = [{[1:100]},{[101:200]}]; % if indices is a double and you want to compare different colors, make each group a separate element in a cell array


% clear datarray
% dat = outdat; 
% latenttype = 'thislatents'; % if empty, then will plot smoothed data
% latentinds = 1:3; % number of latents, must have 2 or 3 elements
% plottype = 'mean' % options: traj, mean, all
% maxtraj = 30; % sets max number of samples, if blank than will do all
% ms = 10; % this is the marker size for scatter plots

    if isequal(class(indices),'double')
        indices = {indices};
    end
    %% extract trials
    trialinds = [{}];
    if isequal('char',class(indices))
        for n = 1:length(indices)
            switch indices(n)
                case 'b'
                    thiscode = 162;
                case 'g'
                    thiscode = 161;
                case 'm'
                    thiscode = 155;
                case 'h'
                    thiscode = 150;
                case 'c'
                    thiscode = [];
            end
            if isempty(thiscode)
                trialinds = 1:dat(1).params.trial.recaltrial(end);
            else
                results = [dat.result];
                trialinds{n} = find(results == thiscode);
            end
        end
    elseif isequal('cell',class(indices))
        trialinds = indices;
    else
        error('Unsupported Index Class')
    end

    colorstouse = distinguishable_colors(length(trialinds));
    h=figure;
    for n = 1:length(trialinds)
        thisinds = trialinds{n}; % get this groups inds
        sample = randperm(length(thisinds)); % sample from them if traj
        if ~isempty(maxtraj)
            if length(thisinds) > maxtraj
                permedinds = sample(1:maxtraj);
            else
                permedinds = sample(1:length(thisinds));
            end
        else
            permedinds = thisinds;
        end
        for m = 1:length(permedinds)
            if isequal('thislatent',latenttype)
                thislatent = dat(permedinds(m)).thislatent;
            else
                thislatent = dat(permedinds(m)).smoothlatent;
            end

            if isequal(plottype,'traj')
                if length(latentinds)==3
                    plot3(thislatent(latentinds(1),:),thislatent(latentinds(2),:),thislatent(latentinds(3),:),'Color',colorstouse(n,:))
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]); zlabel(['Latent ',num2str(latentinds(3))])
                else
                    plot(thislatent(latentinds(1),:),thislatent(latentinds(2),:),'Color',colorstouse(n,:))
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]);                
                end
            elseif isequal(plottype,'mean')
                meanlatents = mean(thislatent,2);
                if length(latentinds)==3
                    scatter3(thislatent(latentinds(1)),thislatent(latentinds(2)),thislatent(latentinds(3)),ms,colorstouse(n,:),'fill')
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]); zlabel(['Latent ',num2str(latentinds(3))])
                else
                    scatter(thislatent(latentinds(1)),thislatent(latentinds(2)),ms,colorstouse(n,:),'fill')
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]);   
                end
            elseif isequal(plottype,'all')
                if length(latentinds)==3
                    scatter3(thislatent(latentinds(1),:),thislatent(latentinds(2),:),thislatent(latentinds(3),:),ms,colorstouse(n,:),'fill')
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]); zlabel(['Latent ',num2str(latentinds(3))])
                else
                    scatter(thislatent(latentinds(1),:),thislatent(latentinds(2),:),ms,colorstouse(n,:),'fill')
                    xlabel(['Latent ',num2str(latentinds(1))]); ylabel(['Latent ',num2str(latentinds(2))]);
                end
            end
            hold on;
        end
    end
    
end

function good_chan_nums = get_good_channels(dat)

    sp_train = [dat.counts];
    good_chan_nums = 1:size(sp_train,1);
    
    % get rid of low fr & high fano
    mean_rates = mean(sp_train,2).*1000;
    low_fr = mean_rates<2;
    
    % get rid of high fano
    minT = min([dat.nBins]);
    count_mat = nan(size(dat(1).counts,1),length(dat));
    for i_trial = 1:length(dat)
        count_mat(:,i_trial) = sum(dat(i_trial).counts(:,1:minT),2);
    end
    
    targ_angs = [dat.angle];
    un_angs = unique(targ_angs);
    fanos = nan(length(good_chan_nums),length(un_angs));
    for i_ang = 1:length(un_angs)
        curr_idx = targ_angs==un_angs(i_ang);
        curr_dat = count_mat(:,curr_idx);
        fanos(:,i_ang) = var(curr_dat,1,2)./mean(curr_dat,2);
    end
    
    avg_fano = mean(fanos,2);
    high_fano = avg_fano>10;
    
    sp_train = sp_train(~low_fr & ~high_fano,:);
    good_chan_nums = good_chan_nums(~low_fr & ~high_fano);
    
    % get rid of cross-talk
    C = normCoincidence(sp_train);
    good_chans = rm_crosstalk_channels(C);
    good_chan_nums = good_chan_nums(good_chans);
    
end


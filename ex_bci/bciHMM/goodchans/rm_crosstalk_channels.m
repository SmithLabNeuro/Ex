function good_chans = rm_crosstalk_channels(coinc_matrix,coinc_thresh)
    
    if nargin<2
        coinc_thresh = 0.2;
    end
    n_chans = size(coinc_matrix,1);

    % figure out good channels
    tmp_coinc_mat = coinc_matrix;
    chan_idx = 1:n_chans;
    [x,y] = find(coinc_matrix>coinc_thresh);
    while ~isempty(x)
        rm_neuron = mode([x; y]);
        tmp_coinc_mat(rm_neuron,:)=[]; tmp_coinc_mat(:,rm_neuron)=[];
        chan_idx(rm_neuron) = [];
        
        [x,y] = find(tmp_coinc_mat>coinc_thresh);
    end
    good_chans = false(n_chans,1);
    good_chans(chan_idx) = true;
    
end


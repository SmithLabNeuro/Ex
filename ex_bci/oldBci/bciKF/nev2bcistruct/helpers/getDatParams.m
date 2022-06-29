function dat = getDatParams(dat)
    for n = 1:length(dat)
    preTrial{n}  = dat(n).text;
    variables = regexp(preTrial{n},';','split');
    for j = 1:length(variables)
        k = strfind(variables{j},'=');
        if k            
            lhs = variables{j}(1:k-1);
            rhs = variables{j}(k+1:end);
            % MATT - should we check here to see if any value has
            % changed and report back if it has?
            try
                eval(['params.' lhs '=[' rhs '];']);
            catch
                eval(['params.' lhs '=''' rhs ''';']);
            end
        end
    end
    dat(n).params.trial = params;
    end
end
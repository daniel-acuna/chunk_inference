function v = learn_variance(chunks, data, gamma, start_pause, nonstart_pause, fit_mean)
% function v = learn_variance(chunks, data, gamma, start_pause, nonstart_pause, fit_mean)
% learn variance similar to learn_pause procedure

ind_chunk_start = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;
n_chunks = size(chunks, 1);
n_seq_len = size(chunks, 2);

goodData = ~isnan(data);   % Indicator for good data 
tmp_var    = zeros(n_chunks, 1); % Sufficient stats for weighted sum-of-squares
tmp_weight = zeros(n_chunks, 1); % Sufficient stats for the summed weights for good data 
for i = 1:n_chunks
    if fit_mean
        % center movement_times
        cmt = bsxfun(@minus, data, ...
            ind_chunk_start(i, :)*start_pause + ...
            (~ind_chunk_start(i, :))*nonstart_pause);
        tmp = bsxfun(@times, cmt.^2, gamma(2:end, i));
    else
        % no centering needed - mean assumed to be zero
        tmp = bsxfun(@times, data.^2, gamma(2:end, i));
    end
    tmp_var(i) = nansum(tmp(:));

    % Keep track of how much data we have for each structure (weighted by
    % Gamma). 
    tmp_weight(i) = sum(sum(bsxfun(@times,goodData,gamma(2:end, i)))); 
end
v = sum(tmp_var)/sum(sum(tmp_weight));
% Comment (JD): this seems to be a unnecessary complication? Since we are
% using all available data, we can just devide by the number of good data points: 
% sum(tmp_var)/sum(sum(goodData));
% Should get indentical results: or do we have a case in which gamma, over
% the different chunking structures does not sum to 1? 
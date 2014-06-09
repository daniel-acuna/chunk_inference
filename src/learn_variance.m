function v = learn_variance(chunks, data, gamma, ...
    start_pause, nonstart_pause, fit_mean)
% learn variance similar to learn_pause procedure

ind_chunk_start = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;
n_chunks = size(chunks, 1);
tmp_var = zeros(n_chunks, 1);
n_seq_len = size(chunks, 2);
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
    tmp_var(i) = sum(tmp(:));
end
v = sum(tmp_var)/(n_seq_len*sum(sum(gamma(2:end, :))));
function [start_pause, nonstart_pause] = learn_pause(chunks, data, gamma)
% learn pause at the start of chunk and pause not at the start
% base on probability of chunks in matrix gamma for the data `data`

% indicator variable for when chunk starts
ind_chunk = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;
sum_mean_pause = 0;
n_mean_pause = 0;
n_chunks = size(chunks, 1);
for i = 1:n_chunks
    tmp_sem = bsxfun(@times, data, ...
        ind_chunk(i, :));
    tmp_sem = bsxfun(@times, tmp_sem, gamma(2:end, i));
    sum_mean_pause = sum_mean_pause + ...
        sum(tmp_sem(:));
    n_mean_pause = n_mean_pause + ...
        sum(nnz(ind_chunk(i, :)) * gamma(2:end, i));
end
start_pause = sum_mean_pause / n_mean_pause;

% indicator variable for when chunk not started
ind_chunk = ~(diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0);
sum_mean_pause = 0;
n_mean_pause = 0;
for i = 1:n_chunks
    tmp_sem = bsxfun(@times, data, ...
        ind_chunk(i, :));
    tmp_sem = bsxfun(@times, tmp_sem, gamma(2:end, i));
    sum_mean_pause = sum_mean_pause + ...
        sum(tmp_sem(:));
    n_mean_pause = n_mean_pause + ...
        sum(nnz(ind_chunk(i, :)) * gamma(2:end, i));
end
nonstart_pause = sum_mean_pause / n_mean_pause;
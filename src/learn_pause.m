function [start_pause, nonstart_pause] = learn_pause(chunks, data, gamma)
% function [start_pause, nonstart_pause] = learn_pause(chunks, data, gamma)
% learn pause at the start of chunk and pause not at the start
% based on probability of chunks in matrix gamma for the data `data`
% indicator variable for when chunk starts
% Ignores NaNs in the data 
% Parameters: 
%   chunks: Possible chunking structures (nChunks x nMovements) 
%   data:   Data to learn from (nTrials x nMovements) 
%   gamma:  Posterior probability of the hidden state (nTrials x nChunks) 
% Returns: 
%   start_pause: Value for beginning of chunk 
%   nonstart_pause: Value for the end of chunk 
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
        nansum(tmp_sem(:));
    % Now check how many observations we have 
    tmp_ind = bsxfun(@times,~isnan(data),ind_chunk(i, :)); 
    tmp_ind = bsxfun(@times,tmp_ind,gamma(2:end, i)); 
    n_mean_pause = n_mean_pause + ...
         sum(tmp_ind(:));
end
start_pause = sum_mean_pause / n_mean_pause;

% indicator variable for when chunk not started
ind_chunk = ~ind_chunk;
sum_mean_pause = 0;
n_mean_pause = 0;
for i = 1:n_chunks
    tmp_sem = bsxfun(@times, data, ...
        ind_chunk(i, :));
    tmp_sem = bsxfun(@times, tmp_sem, gamma(2:end, i));
    sum_mean_pause = sum_mean_pause + ...
        nansum(tmp_sem(:));
    % Now check how many observations we have 
    tmp_ind = bsxfun(@times,~isnan(data),ind_chunk(i, :)); 
    tmp_ind = bsxfun(@times,tmp_ind,gamma(2:end, i)); 
    n_mean_pause = n_mean_pause + ...
         sum(tmp_ind(:));
end
nonstart_pause = sum_mean_pause / n_mean_pause;
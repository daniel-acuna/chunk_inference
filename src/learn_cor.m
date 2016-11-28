function rho = learn_cor(chunks, cor_chunks, data, gamma, ...
    start_pause, nonstart_pause, variance, fit_mean)
% learn correlation

ind_chunk_start = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;

n_numerator = 0;
n_chunks = size(chunks, 1);
n_seq_len = size(chunks, 2);
tmp_cov = zeros(n_chunks, 1);
goodData=~isnan(data); 
for i = 1:n_chunks
    if fit_mean
        % center movement_times
        cmt = bsxfun(@minus, data, ...
            ind_chunk_start(i, :)*start_pause + ...
            (~ind_chunk_start(i, :))*nonstart_pause);
    else
        cmt = mt_seq;
    end
    % Old version: 
    % Sigma = ((bsxfun(@times, cmt, gamma(2:end, i)))'*cmt) .* ...
    %     cor_chunks(:, :, i) .* (1-eye(n_seq_len));
    % n_numerator = n_numerator + ...
    %     sum(gamma(2:end, i))*nnz(cor_chunks(:, :, i) .* (1-eye(n_seq_len)));

    % To ignore the NaN observations: 
    SS  = bsxfun(@times,permute(cmt,[1 2 3]),permute(cmt,[1 3 2]));  % Sums of squares for each individual trial
    wSS = bsxfun(@times,SS,gamma(2:end, i));    % Weighted by the posterior probaility 
    sSS = squeeze(nansum(wSS,1));               % Take sum ignoring the NaNs
    Sigma = sSS .* cor_chunks(:, :, i) .* (1-eye(n_seq_len));
    num   = ((bsxfun(@times, goodData, gamma(2:end, i)))'*goodData) .* ...
        cor_chunks(:, :, i) .* (1-eye(n_seq_len));  
    tmp_cov(i)  = nansum(Sigma(:));
    n_numerator = n_numerator + sum(num(:)); 
end
new_cov = sum(tmp_cov(:))/n_numerator;
rho = new_cov/variance;
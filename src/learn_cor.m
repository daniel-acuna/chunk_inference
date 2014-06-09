function rho = learn_cor(chunks, cor_chunks, data, gamma, ...
    start_pause, nonstart_pause, variance, fit_mean)
% learn correlation

ind_chunk_start = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;

n_numerator = 0;
n_chunks = size(chunks, 1);
n_seq_len = size(chunks, 2);
tmp_cov = zeros(n_chunks, 1);
for i = 1:n_chunks
    if fit_mean
        % center movement_times
        cmt = bsxfun(@minus, data, ...
            ind_chunk_start(i, :)*start_pause + ...
            (~ind_chunk_start(i, :))*nonstart_pause);
    else
        cmt = mt_seq;
    end
    Sigma = ((bsxfun(@times, cmt, gamma(2:end, i)))'*cmt) .* ...
        cor_chunks(:, :, i) .* (1-eye(n_seq_len));
    n_numerator = n_numerator + ...
        sum(gamma(2:end, i))*nnz(cor_chunks(:, :, i) .* (1-eye(n_seq_len)));
    tmp_cov(i) = sum(Sigma(:));
end
new_cov = sum(tmp_cov(:))/n_numerator;
rho = new_cov/variance;
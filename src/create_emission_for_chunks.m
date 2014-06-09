function [mt_cov, p_obs, T, p_obs_un] = ...
    create_emission_for_chunks(rho, rho_er, self_transition, mt_seq, ...
    er_seq, varargin)
% Create emission structure for chunks.
% Create probability of observation and movement time covariances 
% according to chunk.
%
% Parameters:
%    rho        : correlation within chunks
%    rho_er     : correlation of errors
%    self_t     : self transition probability
%    mt_seq     : detrended movement times observed
%    er_seq     : detrended error sequence
%    mean_pause : pause at beginning of chunk
%    
% Returns:
%    mt_cov: movement time covariance for all chunks
%    p_obs : probability of observation
%    T     : transition probabilities across chunks

p = inputParser;
p.KeepUnmatched = true;
addParamValue(p, 'v', var(mt_seq(:)), @isnumeric);
addParamValue(p, 'v_er', var(er_seq(:)), @isnumeric);
% addParamValue(p, 'model', 'mt_er', @(x)ismember(x, {'mt', 'er', 'mt_er'}));
addParamValue(p, 'fit_er', true, @islogical);
addParamValue(p, 'fit_rt', true, @islogical);
addParamValue(p, 'mean_pause', 0, @isnumeric);
addParamValue(p, 'mean_inchunk', 0, @isnumeric);
addParamValue(p, 'mean_pause_er', 0, @isnumeric);
addParamValue(p, 'mean_inchunk_er', 0, @isnumeric);
addParamValue(p, 'chunks', [], @ismatrix);
parse(p, varargin{:});
v = p.Results.v;
v_er = p.Results.v_er;
fit_er = p.Results.fit_er;
fit_rt = p.Results.fit_rt;
mean_pause = p.Results.mean_pause;
mean_inchunk = p.Results.mean_inchunk;
mean_pause_er = p.Results.mean_pause_er;
mean_inchunk_er = p.Results.mean_inchunk_er;
chunks = p.Results.chunks;

if ~fit_er && ~fit_rt
    error('Either fit_er or fit_mt must be true');
end
% Find best parameters for chunking

if isempty(chunks)
    chunks = create_chunks3();    
end
cor_chunks = to_corr_chunks(chunks);
n_chunks = size(chunks, 1);

[chunk_means_mt, mt_cov, chunk_means_er, er_cov] = ...
    create_chunk_means_covs(chunks, cor_chunks, ...
        mean_pause, mean_inchunk, v, rho, ...
        mean_pause_er, mean_inchunk_er, v_er, rho_er);

% Compute probability of observation

n_time = size(mt_seq, 1);
p_obs = zeros(n_time, n_chunks);

for i = 1:n_chunks
    if fit_rt && fit_er    
        p_obs(:, i) = ...
            gaussLogprob(chunk_means_mt(i, :)', ...
                    mt_cov(:, :, i), mt_seq) + ...
            gaussLogprob(chunk_means_er(i, :)', ...
                    er_cov(:, :, i), er_seq);
    elseif fit_rt
        p_obs(:, i) = ...
            gaussLogprob(chunk_means_mt(i, :)', ...
            mt_cov(:, :, i), mt_seq);
    else
        p_obs(:, i) = gaussLogprob(chunk_means_er(i, :)', ...
            er_cov(:, :, i), er_seq);
    end
end
% p_obs = exp(normalizeLogspace(p_obs));
p_obs = exp(p_obs);
% Unnormalized
p_obs_un = p_obs;

% Normalize
p_obs = bsxfun(@rdivide, p_obs, sum(p_obs, 2));

% Transition probability probability of self-transition
T = ((ones(n_chunks, n_chunks)-eye(n_chunks))*...
    (1-self_transition))/(n_chunks-1) + ...
    eye(n_chunks)*self_transition;

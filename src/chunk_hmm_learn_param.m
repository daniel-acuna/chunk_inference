function [rho, self_t, log_like, fm, T, rho_er, v, v_er, ...
    initial_dist, mean_pause, mean_inchunk, ...
    mean_pause_er, mean_inchunk_er, ...
    chunks, cor_chunks, gamma] = ...
    chunk_hmm_learn_param(mt_seq, er_seq, varargin)
% Finds best parameters for chunks using Baum-Welch algorithm
% input parameters:
%  mt_seq: reaction times n x k (n: number of trials, k: sequence length)
%  er_seq: errors n x k
%  rho   : initial correlation of reaction times within chunks
%  rho_er: initial correlation of errors within chunks
%  self_t: initial self transition probabilities
%  model : specify whether to use movement times (mt), errors (er), or 
%          both (mt_er) for inferring chunking
%  fit_pause: specify whether a pause is expected at thebeginning of a 
%          chunk 
%  verbose:whether to show the fitting process
%  mean_pause:specify the initial mean pause expected

% output parameters
% rho: best correlation parameter for reaction time
% self_t: best self-transition parameter
% log_like: log-likelihood on training data
% fm: forward messages for hmm learning
% T: transition matrix
% rho_er : correlation parmater for errors within a chunk
% v: variance in reaction time
% v_er: variance in errors
% initial_dist: initial distribution of chunks
% mean_pause: mean reaction time pause at beginning of chunk

% some parameters can be given
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'mt_seq', @ismatrix);
addRequired(p, 'er_seq', @(x)(all(size(mt_seq) == size(er_seq))));
addParamValue(p, 'fit_rt', false, @islogical);
addParamValue(p, 'fit_er', false, @islogical);
addParamValue(p, 'fit_rt_rt', false, @islogical);
addParamValue(p, 'fit_er_er', false, @islogical);
addParamValue(p, 'fit_T', true, @islogical);
addParamValue(p, 'fit_rho', false, @islogical);
addParamValue(p, 'fit_rho_er', false, @islogical);
addParamValue(p, 'diagonal_T', false, @islogical);
addParamValue(p, 'rho', 0.1, @isscalar);
addParamValue(p, 'rho_er', 0.1, @isscalar);
addParamValue(p, 'self_t', 0.6, @isscalar);
addParamValue(p, 'verbose', false, @islogical);
addParamValue(p, 'mean_pause', 0.1, @isnumeric);
addParamValue(p, 'mean_inchunk', 0, @isnumeric);
addParamValue(p, 'mean_pause_er', 0.1, @isnumeric);
addParamValue(p, 'mean_inchunk_er', 0, @isnumeric);
addParamValue(p, 'chunks', []);
addParamValue(p, 'cor_chunks', []);
parse(p, mt_seq, er_seq, varargin{:});
rho = p.Results.rho;
rho_er = p.Results.rho_er;
self_t = p.Results.self_t;
fit_rt = p.Results.fit_rt;
fit_er = p.Results.fit_er;
fit_rt_rt = p.Results.fit_rt_rt;
fit_er_er = p.Results.fit_er_er;
fit_T = p.Results.fit_T;
fit_rho = p.Results.fit_rho;
fit_rho_er = p.Results.fit_rho_er;
diagonal_T = p.Results.diagonal_T;
verbose = p.Results.verbose;
mean_pause = p.Results.mean_pause;
mean_inchunk = p.Results.mean_inchunk;
mean_pause_er = p.Results.mean_pause_er;
mean_inchunk_er = p.Results.mean_inchunk_er;
chunks = p.Results.chunks;

tic;

if ~fit_rt
    mean_pause = 0;
    mean_inchunk = 0;
end

if ~fit_er
    mean_pause_er = 0;
    mean_inchunk_er = 0;
end

if ~fit_rt_rt
    rho = 0;
end

if ~fit_er_er
    rho_er = 0;
end

log_like = nan;

if isempty(chunks)
    chunks = create_chunks_nospace('n_seqlen', size(mt_seq, 2));
end

cor_chunks = to_corr_chunks(chunks);


n_chunks = size(chunks, 1);
T = ((ones(n_chunks, n_chunks)-eye(n_chunks))*...
    (1-self_t))/(n_chunks-1) + ...
    eye(n_chunks)*self_t;

not_exit = true;
v = var(mt_seq(:));
v_er = var(er_seq(:));

% initial distribution
initial_dist = ones(1, n_chunks)/n_chunks;
n_iteration = 0;
% alpha = 1;
% beta = 1/(n_chunks-1);
while not_exit    
    current_log_like = sum(log_like);
    % COMPUTATION OF EXPECTATION
    [~, ~, T0, p_obs_un] = ...
        create_emission_for_chunks(rho, rho_er, self_t, mt_seq, ...
            er_seq, 'v', v, 'v_er', v_er, ...
            'fit_er', fit_er || fit_er_er, ...
            'fit_rt', fit_rt || fit_rt_rt, ...
            'mean_pause', mean_pause, ...
            'mean_inchunk', mean_inchunk, ...
            'mean_pause_er', mean_pause_er, ...
            'mean_inchunk_er', mean_inchunk_er, ...
            'chunks', chunks, ...
            'cor_chunks', cor_chunks);
    if ~fit_T
        T = T0;
    end
    
%     % Small prior
%     prior = alpha*eye(size(T)) + beta*(ones(size(T)) - eye(size(T)));
%     T = T + prior;
%     % Normalize
%     T = bsxfun(@rdivide, T, sum(T, 2));
%     disp(self_t);

    [fm, ~, gamma, log_like, marg_epsilon] = ...
        hmm_inference(p_obs_un, T, 'initial_dist', initial_dist);
    if verbose
        n_iteration = n_iteration + 1;
        
        fprintf('Iteration: %.0f (elapsed time: %.1f s)\n', n_iteration, ...
            toc);
        fprintf(...
            ['\t rho             = %f\n' ...
             '\t rho_er          = %f\n' ...
             '\t self_t          = %f\n' ...
             '\t v               = %f\n' ...
             '\t v_er            = %f\n' ...
             '\t mean_pause      = %f\n' ...
             '\t mean_inchunk    = %f\n' ...
             '\t mean_pause_er   = %f\n' ...
             '\t mean_inchunk_er = %f\n' ...
             '\t log_like        = %f\n'], ...
            rho, rho_er, self_t, v, v_er, ...
            mean_pause, mean_inchunk, ...
            mean_pause_er, mean_inchunk_er, sum(log_like));
        %disp(sum(log_like));
        %disp(self_t);
%         plot(gamma);
%         drawnow;
    end
    delta_fval = abs(sum(log_like) - current_log_like);
    avg_fval = ((abs(sum(log_like)) + abs(current_log_like))+eps)/2;        
    if ~isnan(current_log_like) && ...
            ((delta_fval/avg_fval < 1e-8) || ...
            (sum(log_like) - current_log_like < - 2*eps))
        
        % if log-likelihood went down, then return previous result
        % save previous results
        if (sum(log_like) - current_log_like < - 2*eps)
            rho = rho_best;
            self_t = self_t_best;
            T = T_best;
            rho_er = rho_er_best;
            v = v_best;
            v_er = v_er_best;
            initial_dist = initial_dist_best;
            mean_pause = mean_pause_best;
            mean_inchunk = mean_inchunk_best;
            mean_pause_er = mean_pause_er_best;
            mean_inchunk_er = mean_inchunk_er_best;
            gamma = gamma_best;
            log_like = log_like_best;
        end
        not_exit = false;
    else
        % save previous results
        rho_best = rho;
        self_t_best = self_t;
        T_best = T;
        rho_er_best = rho_er;
        v_best = v;
        v_er_best = v_er;
        initial_dist_best = initial_dist;
        mean_pause_best = mean_pause;
        mean_inchunk_best = mean_inchunk;
        mean_pause_er_best = mean_pause_er;
        mean_inchunk_er_best = mean_inchunk_er;
        gamma_best = gamma;
        log_like_best = log_like;
        
        % MAXIMIZATION OF EXPECTATION
        % Restimate transitions
        % maximization
        if fit_T
            T = marg_epsilon;

            self_t = trace(gamma(1:end-1, :)' * ...
                gamma(2:end, :))/sum(sum(gamma(1:end-1, :)' * ...
                    gamma(2:end, :)));
        end
        if fit_T && diagonal_T
            % Transition probability probability of self-transition
            T = ((ones(n_chunks, n_chunks)-eye(n_chunks))*...
                (1-self_t))/(n_chunks-1) + ...
                eye(n_chunks)*self_t;
        end
        
        initial_dist = gamma(1, :);        
        
        % reestimate variance and covariance
        if (fit_rt || fit_rt_rt)
            % find mean_pause
            if fit_rt
                [mean_pause, mean_inchunk] = ...
                    learn_pause(chunks, mt_seq, gamma);
            end
            
            v = learn_variance(chunks, mt_seq, gamma, ...
                mean_pause, mean_inchunk, fit_rt);

            % reestimate covariance for movement time
            if fit_rho
                rho = learn_cor(chunks, cor_chunks, mt_seq, gamma, ...
                    mean_pause, mean_inchunk, v, fit_rt);
            end
        end
        
        if (fit_er || fit_er_er)
            % find mean_pause
            if fit_er
                [mean_pause_er, mean_inchunk_er] = ...
                        learn_pause(chunks, er_seq, gamma);
            end
            
            v_er = learn_variance(chunks, er_seq, gamma, ...
                mean_pause_er, mean_inchunk_er, fit_er);

            if fit_rho_er
                rho_er = learn_cor(chunks, cor_chunks, er_seq, gamma, ...
                    mean_pause_er, mean_inchunk_er, v_er, fit_er);
            end
        end
    end
end


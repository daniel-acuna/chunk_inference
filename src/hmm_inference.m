function [fm, bm, gamma, log_like, marg_epsilon] = hmm_inference(p_obs, T, ...
    varargin)
% Infer forward and backward probability of chunks
% and the smoothed estimates
% Inputs:
%   p_obs : probability of observations
%   T     : transition probability from i to j
%   fm0   : initial distribution on chunks (optional)
% Outputs:
%   fm     : forward messages
%   bm     : backward messages
%   gamma  : smoothed estimates
%   epsilon: 

p = inputParser;
addRequired(p, 'p_obs', @ismatrix);
addRequired(p, 'T', @ismatrix);
addParamValue(p, 'initial_dist', [], @ismatrix);
parse(p, p_obs, T, varargin{:})
initial_dist = p.Results.initial_dist;

[n_time, n_chunks] = size(p_obs);

if (isempty(initial_dist))
    initial_dist = ones(1, n_chunks)/n_chunks;
end   

% normalizing factors
Z = ones(n_time+1, 1);

% forward messages
fm = zeros(n_time + 1, n_chunks);
% flat prior
fm(1, :) = initial_dist;


for t = 2:(n_time+1)
    % Compute probability of observed reaction times and errors
    % given the chunks
    % Forward message
    fm(t, :) = p_obs(t-1, :)'.*(T'*fm(t-1, :)');
    Z(t) = sum(fm(t, :));
    fm(t, :) = fm(t, :)/Z(t);
end

% Compute backward messages
% disp('Backward pass...');
% Backward messages
bm = zeros(n_time+1, n_chunks);
bm(end, :) = 1;
for t = (n_time+1):-1:2
    bm(t-1, :) = T*(p_obs(t-1, :)'.*bm(t, :)');
    bm(t-1, :) = bm(t-1, :)/Z(t);
end

% compute smoothed estimations
gamma = fm .* bm;

% compute expected probabilities
% epsilon = zeros(n_chunks, n_chunks, n_time-1);
%mepsilon = zeros(n_time-1, n_chunks, n_chunks);
% normalize p_obs
% p_obs_norm = bsxfun(@rdivide, p_obs, sum(p_obs, 2));
% % Marginalized epsilon across time
% % marg_epsilon = (((gamma(1:(n_time - 1), :) ./ bm(1:(n_time - 1), :)) * T)' * ...
% %     (p_obs_norm(2:n_time, :) .* bm(2:n_time, :)))';
% marg_epsilon = ((gamma(1:end-1, :) ./ bm(1:end-1, :)) * T)' * ...
%     (p_obs_norm .* bm(2:end, :));
% for i = 1:n_chunks
%     disp(i);
%     tmp = bsxfun(@times, ...
%         gamma(1:(n_time-1), :) ./ bm(1:(n_time - 1), :), T(:, i)');
%     epsilon(:, :, i) = bsxfun(@times, ...
%             tmp, ...
%             p_obs(2:n_time, i).*bm(2:n_time, i));
% end

% for t = 1:(n_time - 1)
%     epsilon2(t, :, :) = epsilon2
%     bsxfun(@times, epsilon2, sum(epsilon2, 3));
% end

% epsilon2 = shiftdim(epsilon2, 2);
% 
if n_chunks < 200
%     epsilon = zeros(n_chunks, n_chunks, n_time);
%     for i = 1:n_chunks
%         for j = 1:n_chunks
% %             epsilon(i, j, :) = ((gamma(1:(n_time-1), i).*T(i, j).*...
% %                 p_obs(2:n_time, j).*bm(2:n_time, j)) ./ ...
% %                 bm(1:(n_time - 1), i));
%             epsilon(i, j, :) = ((gamma(1:end-1, i).*T(i, j).*...
%                 p_obs(:, j).*bm(2:end, j)) ./ ...
%                 bm(1:end-1, i));
%         end
%     end
% 
%     % Normalize
%     for t = 1:(n_time-1)
%         epsilon(:, :, t) = epsilon(:, :, t) / sum(sum(epsilon(:, :, t)));
%     end
%     marg_epsilon2 = sum(epsilon, 3);    
end

% How to get the same value as marg_epsilon2
p_obs_norm = bsxfun(@rdivide, p_obs, sum(p_obs, 2));
% Marginalized epsilon across time
% marg_epsilon = (((gamma(1:(n_time - 1), :) ./ bm(1:(n_time - 1), :)) * T)' * ...
%     (p_obs_norm(2:n_time, :) .* bm(2:n_time, :)))';
%bm_norm = (p_obs_norm .* bm(2:end, :));

% bm_norm = bsxfun(@rdivide, bm(2:end, :), sum(bm(2:end, :), 2));

marg_epsilon = ((gamma(1:end - 1, :) ./ bm(1:end - 1, :))' * ...
    (p_obs_norm .* bm(2:end, :))) .* T;
marg_epsilon = bsxfun(@rdivide, marg_epsilon, sum(gamma)');
% marg_epsilon = ((gamma(1:end - 1, :) ./ bm(1:end - 1, :))' * ...
%     (p_obs_norm .* bm_norm)) .* T;
% marg_epsilon = ((fm(1:end - 1, :)' * (p_obs_norm .* bm(1:end-1, :))) .* T);
marg_epsilon = bsxfun(@rdivide, marg_epsilon, sum(marg_epsilon, 2));
% marg_epsilon2 = bsxfun(@rdivide, marg_epsilon2, sum(marg_epsilon2, 2));
% 
% T0 = (gamma(1:end-1, :) ./ bm(1:end-1, :)) * T;
% T0 = bsxfun(@rdivide, T0, sum(T0, 2));
% marg_epsilon = T0' * ...
%     (p_obs_norm .* bm(2:end, :));


log_like = (log(Z(2:end)));
% log_like = sum(log(sum(gamma(2:end, :) .* p_obs, 2)));

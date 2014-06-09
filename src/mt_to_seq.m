function [mt_seq, er_seq] = mt_to_seq(mt, mt_residuals, ...
    er_residuals, varargin)
% Return matrix of sequence_press vs keypress
% for movement time and error

p = inputParser;
addParamValue(p, 'remove_nan', true, @islogical);
parse(p, varargin{:});
remove_nan = p.Results.remove_nan;


n_seq = length(unique(mt.sequence_trial));
n_key = length(unique(mt.sequence_press));

mt_seq = nan(n_seq, n_key);
er_seq = nan(n_seq, n_key);
movement_time = double(mt_residuals);
%errors = double(mt(:, 'error'));
errors = double(er_residuals);
trial = double(mt.sequence_trial);
keypress = double(mt.sequence_press);
all_trials = double(unique(mt.sequence_trial));
all_keys = double(unique(mt.sequence_press));
for i = 1:n_seq
    %day_id(i) = days(find(trial == all_trials(i), 1));
    for j = 1:n_key
        idx = find(trial == all_trials(i) & keypress == all_keys(j), 1);
        if ~isempty(idx)
            mt_seq(i, j) = movement_time(idx);
            er_seq(i, j) = errors(idx);
        end
    end
end

if remove_nan
    % Get ids that don't have NA
    idx = ~((any(isnan(mt_seq), 2)) | (any(isnan(er_seq), 2)));

    mt_seq = mt_seq(idx, :);
    er_seq = er_seq(idx, :);
end
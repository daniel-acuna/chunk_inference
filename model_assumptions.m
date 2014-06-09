d = load('full_data/mt.mat');
mt = d.mt;
mt = mt(mt.sequence_trial <= 1800 & mt.sequence_press <= 10 & ...
    mt.movement_time >= 0.05 & mt.movement_time <= 2, :); %#ok<NODEF>
%%

res = {};
res_model = {};
for i = unique(mt.subject_id)'
    for j = unique(mt.sequence_id)'
        nlmf = NonLinearModel.fit(...
            mt(mt.subject_id == i & mt.sequence_id == j, :), ...
            exponential_model, initial_values, 'Options', opts);
        res{end+1} = [nlmf.Residuals(:, 'Raw') ...
            mt(mt.subject_id == i & mt.sequence_id == j, :)];
        res_model{end+1} = ...
            LinearModel.fit(res{end}, ...
            'Raw ~ sequence_press + sequence_trial');
    end
end

%% testing hypothesis that RRT is independent of sequence trial and
% sequence press
[h, p] = ttest(cellfun(@(x)(corr(x.sequence_press, x.Raw)), res))

[h, p] = ttest(cellfun(@(x)(corr(x.sequence_trial, x.Raw)), res))

%% effect size
mean(cell2mat(cellfun(@(x)(double(x.Coefficients(:, 'Estimate'))), ...
    res_model, 'UniformOutput', false)), 2).*[1 1000 1]'

% Range of values
mean(cell2mat(cellfun(@(x)(quantile(x.Raw, [0.025 1-0.025])'), res, 'UniformOutput', false)), 2)

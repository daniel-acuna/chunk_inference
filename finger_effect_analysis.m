% removing effect of finger on human data
addpath('src/');
load('full_data/mt_with_fingers.mat');
%%
mt = mt(mt.sequence_trial <= 1800 & mt.sequence_press <= 10 & ...
    mt.movement_time >= 0.05 & mt.movement_time <= 2, :); %#ok<NODEF>

%% analyze all subjects
for subject = head(unique(mt.subject_id), 1)'
    for sequence_id = 1:2
        mt_ss = mt(mt.subject_id == subject & ...
            mt.sequence_id == sequence_id, :);
    end
end

%%
%mt_ss.target = nominal(mt_ss.target);
%mt_ss.target2 = mt_ss.target == 2;
control_for_finger = false;
% Try model that controlls for finger differences
switch control_for_finger
    case true
        exponential_model = ['movement_time ~ a0' ...
            '+ a1*exp((b1/100)*(sequence_trial-1) + ' ...
            'b2/10*(sequence_trial-1)*(sequence_press-1)) +' ...
            'f2 * (target == 2) + ' ...
            'f3 * (target == 3) + ' ...
            'f4 * (target == 4) + ' ...
            'f5 * (target == 5)'];
        
        initial_values = [0.18 0.38 -0.17, -0.1, 0, 0, 0, 0];
        
    case false
        exponential_model = ['movement_time ~ a0' ...
            '+ a1*exp((b1/100)*(sequence_trial-1) + ' ...
            'b2/10*(sequence_trial-1)*(sequence_press-1))'];
        initial_values = [0.18 0.38 -0.17, -0.1];
end


opts = statset('Display','off','TolFun',1e-5, ...
    'MaxIter', 100);

% Fit
nlmf = NonLinearModel.fit(mt_ss, ...
    exponential_model, initial_values, 'Options', opts);

%% Plot detrended response times
figure(2);
clf;
plot(mt_ss.sequence_trial, nlmf.Residuals(:, 'Raw'), '.');
hold on;
plot(mt_ss.sequence_trial, ...
    smooth(double(nlmf.Residuals(:, 'Raw')), 500), 'r-');
xlabel('Trial');
ylabel('Residual response time');
%%

% Create a response time and error matrix of trials vs element 
[rt_seq, er_seq] = mt_to_seq(mt_ss, ...
                             double(nlmf.Residuals(:, 'Raw')), ...
                             mt_ss.error);
%%
% Plot results
figure(3);
clf;
subplot(2, 1, 1);
imagesc(rt_seq', [-0.25 0.25]);
colormap('cool');
colorbar;
xlabel('Trial');
ylabel('Element');
title('Response times');
subplot(2, 1, 2);
imagesc(er_seq', [0 1]);
colorbar;
xlabel('Trial');
ylabel('Element');
title('Errors');

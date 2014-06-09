%% Load raw data
addpath('src/');
load('example_data/dsp_example.mat');

%% Preliminary visualizations
figure(1);
clf;
subplot(2, 1, 1);
plot(rt_er_data.sequence_trial, rt_er_data.movement_time, '.');
hold on;
h = plot(rt_er_data.sequence_trial, smooth(rt_er_data.movement_time, 500), 'r-');
legend(h, 'smoothed reponse time');
xlabel('Trial');
ylabel('Seconds');
title('Reponse time');
subplot(2, 1, 2);
plot(rt_er_data.sequence_trial, rt_er_data.error, '.');
hold on;
h = plot(rt_er_data.sequence_trial, smooth(rt_er_data.error, 500), 'r-');
legend(h, 'smoothed error rate');
xlabel('Trial');
ylabel('Error');
title('Errors');

%% Data detrending
% Model
exponential_model = ['movement_time ~ a0' ...
    '+ a1*exp((b1/100)*(sequence_trial-1) + ' ...
    'b2/10*(sequence_trial-1)*(sequence_press-1))'];
initial_values = [0.18 0.38 -0.17, -0.1];

opts = statset('Display','off','TolFun',1e-5, ...
    'MaxIter', 100);

% Fit
nlmf = NonLinearModel.fit(rt_er_data, ...
    exponential_model, initial_values, 'Options', opts);

% Plot detrended response times
figure(2);
clf;
plot(rt_er_data.sequence_trial, nlmf.Residuals(:, 'Raw'), '.');
hold on;
plot(rt_er_data.sequence_trial, ...
    smooth(double(nlmf.Residuals(:, 'Raw')), 500), 'r-');
xlabel('Trial');
ylabel('Residual response time');
                         
%% Response time and error matrices

% Create a response time and error matrix of trials vs element 
[rt_seq, er_seq] = mt_to_seq(rt_er_data, ...
                             double(nlmf.Residuals(:, 'Raw')), ...
                             rt_er_data.error);

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

%% Create space of chunk structures
chunk_structures = create_chunks_nospace('n_seqlen', size(rt_seq, 2));
figure(4);
clf;
% Plot space of chunks
imagesc(chunk_structures');
xlabel('Chunk structure index');
ylabel('Element');
colormap('jet');
title('All possible chunk structures');

%% Find which chunk structure is present at each trial
[rho, self_t, log_like, fm, T, rho_er, v, v_er, ...
    initial_dist, mean_pause, mean_inchunk, ...
    mean_pause_er, mean_inchunk_er, ...
    chunks, cor_chunks, gamma] = ...
    chunk_hmm_learn_param(rt_seq, er_seq, 'verbose', true, ...
    'fit_rt', true, 'fit_rt_rt', true, 'fit_er', true, 'fit_er_er', true, ...
    'fit_T', true, 'fit_rho', true, 'fit_rho_er', true, ...
    'chunks', chunk_structures);

% compute mean and covariance of each chunk structure
[chunk_means_rt, rt_cov, chunk_means_er, er_cov] = ...
    create_chunk_means_covs(chunks, cor_chunks, ...
        mean_pause, mean_inchunk, v, rho, ...
        mean_pause_er, mean_inchunk_er, v_er, rho_er);

%% Plot results of algorithm

% Expected chunking structure
figure(5);
clf;
subplot(2, 1, 1);
imagesc((gamma * chunks)');
colormap('jet');
xlabel('Trial');
ylabel('Element');
title('Expected chunking structure');
subplot(2, 1, 2);

n_chunks = apply(@(x)(length(unique(x))), chunks);

% Expected number of chunk per trial (with some smoothing)
plot(smooth(gamma * n_chunks, 100), '-')
xlabel('Trial');
ylabel('Number of chunks');
title('Expected number of chunks per trial');
axis tight;

% Meean response times and errors fitted by model
expected_rt = gamma * chunk_means_rt;
expected_er = gamma * chunk_means_er;
figure(6);
clf;
subplot(2, 1, 1);
imagesc(expected_rt', [-0.25 0.25]);
colormap('cool');
colorbar;
xlabel('Trial');
ylabel('Element');
title('Expected response times');
subplot(2, 1, 2);
imagesc(expected_rt', [0 0.2]);
colorbar;
xlabel('Trial');
ylabel('Element');
title('Expected error rate');






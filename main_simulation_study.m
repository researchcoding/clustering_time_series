% %% Generate weakly stationary paths
% clear; clc;
% rng(2018)
% % define simulation parameters
% n_clusters = 5;
% obs_num_per_cluster = 10;
% H = -0.4:0.2:0.4;
% total_time_steps = 100;
% obs_num_per_step = 3;
% total_num_observations = total_time_steps * obs_num_per_step;
% total_num_paths = 100;
% 
% % simulation weakly stationary stochastic processes
% [obs_chain, cluster_ind] = sim_wssp_paths(n_clusters, obs_num_per_cluster, H, ...
%                          total_num_observations, total_num_paths);
% save('mBm.mat')
                     
% %% Offline dataset experiments
test_time_steps = 15; 
test_num_sims = 100; 
miscls_rate_offline_algo1 = zeros(test_time_steps, test_num_sims);
miscls_rate_offline_algo2 = zeros(test_time_steps, test_num_sims);
avg_miscls_rate_offline_algo1 = zeros(test_time_steps,1);
avg_miscls_rate_offline_algo2 = zeros(test_time_steps,1);

take_log = 1;
weight_squared = 1;

for t = 1:test_time_steps
    parfor sim = 1:test_num_sims
    % parfor sim = 1:test_num_sims;  % parallel computing if necessary
        
        % scale the obsersed times series to be mean 0
        obs = obs_chain(1:(t * obs_num_per_step), :, sim)';
        obs = scale_mean(obs, 0);
        
        % full matrix is observed under offline dataset
        obs_idx = ones(size(obs));
        
        % clustering the observed time series
        [I_chain_algo1, dm] = unsup_wssp_offline_algo(obs, obs_idx, n_clusters, take_log, weight_squared);
        [I_chain_algo2, ~] = unsup_wssp_online_algo(obs, obs_idx, n_clusters, take_log, weight_squared, dm);
       
        % calculate misclassification rate
        miscls_rate_offline_algo1(t,sim) = misclassify_rate(I_chain_algo1, cluster_ind);
        miscls_rate_offline_algo2(t,sim) = misclassify_rate(I_chain_algo2, cluster_ind);

        fprintf('Offline simluation iter %i for time step %i. \n', sim, t)
    end
    avg_miscls_rate_offline_algo1(t) = mean(miscls_rate_offline_algo1(t,:));
    avg_miscls_rate_offline_algo2(t) = mean(miscls_rate_offline_algo2(t,:));
    
    if mod(t,5) == 0
        save('offline_results_with_log.mat')
    end
end

% plot of clsutering results
x = 1:test_time_steps;
figure
plot(x, avg_miscls_rate_offline_algo1(1:test_time_steps), 'b', 'LineWidth', 2)
hold on
plot(x, avg_miscls_rate_offline_algo2(1:test_time_steps), '-.r', 'LineWidth', 2)
hold off
title('Offline Dataset with Covariance Distance Clustering')
xlabel('time step')
ylabel('misclassification rate')
legend('Algorithm 1', 'Algorithm 2')        
save('offline_results_no_log.mat')


%% Online dataset experiments
test_num_sims = 100; 
miscls_rate_online_algo1 = zeros(test_time_steps, test_num_sims);
miscls_rate_online_algo2 = zeros(test_time_steps, test_num_sims);
avg_miscls_rate_online_algo1 = zeros(test_time_steps,1);
avg_miscls_rate_online_algo2 = zeros(test_time_steps,1);

for t = 1:test_time_steps
    parfor sim = 1:test_num_sims
    % parfor sim = 1:test_num_sims;  % parallel computing if necessary
        
        % scale the obsersed times series to be mean 0
        obs = obs_chain(1:(t * obs_num_per_step), :, sim)';
        
        % full matrix is observed under offline dataset
        obs_idx = zeros(size(obs));
        for i = 1:obs_num_per_cluster:(n_clusters * obs_num_per_cluster)
            obs_idx(i:(i+4), :) = 1;
            for j = 1:(obs_num_per_cluster - 5)
                if t > j * 10
                    obs_idx(5+j, (j * 10 * obs_num_per_step + 1):end) = 1;
                end
            end
        end
        
        keep_idx = find(sum(obs_idx,2) ~= 0);
        obs = obs(keep_idx, :);
        obs_idx = obs_idx(keep_idx, :);
        cluster_ind_online = cluster_ind(keep_idx);
        
        % clustering the observed time series
        [I_chain_algo1, dm] = unsup_wssp_offline_algo(obs, obs_idx, n_clusters, take_log, weight_squared);
        [I_chain_algo2, ~] = unsup_wssp_online_algo(obs, obs_idx, n_clusters, take_log, weight_squared, dm);
       
        % calculate misclassification rate
        miscls_rate_online_algo1(t,sim) = misclassify_rate(I_chain_algo1, cluster_ind_online);
        miscls_rate_online_algo2(t,sim) = misclassify_rate(I_chain_algo2, cluster_ind_online);

        fprintf('Online simluation iter %i for time step %i. \n', sim, t)
    end
    avg_miscls_rate_online_algo1(t) = mean(miscls_rate_online_algo1(t,:));
    avg_miscls_rate_online_algo2(t) = mean(miscls_rate_online_algo2(t,:));
    
    if mod(t,5) == 0
        save('online_results_with_log.mat')
    end
end

% plot of clsutering results
x = 1:test_time_steps;
figure
plot(x, avg_miscls_rate_online_algo1, 'b', 'LineWidth', 2)
hold on
plot(x, avg_miscls_rate_online_algo2, '-.r', 'LineWidth', 2)
hold off
title('Online Dataset with Covariance Distance Clustering')
xlabel('time step')
ylabel('misclassification rate')
legend('Algorithm 1', 'Algorithm 2')     
% save('online_results_no_log.mat')
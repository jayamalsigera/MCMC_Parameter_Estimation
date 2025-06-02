% === MCMC Parameter Estimation for Rotary Inverted Pendulum ===
% using the Metropolis-Hastings algorithm.

clear; clc;

%% === Add Path and Load Simulated Data ===
addpath('D:\MATLAB\MCMC_Parameter_Estimation\Inverted_Pendulum_Model');
load('rip_sim_data.mat');  % loads: t, x

%% === Initial Setup ===
param_names = {'Br', 'Bp', 'kt', 'km', 'eta_m', 'eta_g'};
x0 = [0; 0; pi - 0.1; 0];  % initial state
Vm = @(t) chirp(t, 0.1, 100, 5, 'linear');  % input voltage

%% === Gaussian Priors (mean and std from similar motors) ===
mu_prior = [0.002, 0.002, 0.005, 0.005, 0.7, 0.8];     % based on prior knowledge
sigma_prior = [0.001, 0.001, 0.002, 0.002, 0.3, 0.3];  % wide prior to reflect uncertainty

%% === MCMC Settings ===
num_iters = 1000;
param_dim = length(mu_prior);
burn_in = floor(num_iters / 2);  % burn-in period

% Proposal step sizes
step_size = 0.05 * sigma_prior;

% Bounds for parameters
param_bounds = [1e-4, 1e-2;
                1e-4, 1e-2;
                1e-4, 1e-2;
                1e-4, 1e-2;
                0.01, 1.0;
                0.01, 1.0];

%% === Initialization ===
params_current = mu_prior;
x_ref = x;  % reference data
chain = zeros(num_iters, param_dim);
errors = zeros(num_iters, 1);

% Initial log posterior
x_sim_current = simulate_system(params_current, x0, t, Vm);
error_current = x_ref(:,3) - x_sim_current(:,3);
noise_std = 0.05;
log_like_current = -0.5 * sum((error_current / noise_std).^2);
log_prior_current = -0.5 * sum(((params_current - mu_prior)./sigma_prior).^2);
log_post_current = log_like_current + log_prior_current;

%% === MCMC Sampling Loop ===
tic;
for i = 1:num_iters
    fprintf('[%s] MCMC Iteration %d / %d\n', datestr(now, 'HH:MM:SS'), i, num_iters);

    % Propose new parameters
    proposal = params_current + step_size .* randn(1, param_dim);

    % Reject if out of bounds
    if any(proposal < param_bounds(:,1)') || any(proposal > param_bounds(:,2)')
        chain(i,:) = params_current;
        errors(i) = error_current;
        continue;
    end

    % Compute log-posterior of proposal
    error_prop = compute_error(proposal, x0, t, Vm, x_ref);
    log_like_prop = -0.5 * (error_prop / noise_std)^2;
    log_prior_prop = -0.5 * sum(((proposal - mu_prior)./sigma_prior).^2);
    log_post_prop = log_like_prop + log_prior_prop;

    % Accept/Reject step
    if log(rand) < (log_post_prop - log_post_current)
        params_current = proposal;
        log_post_current = log_post_prop;
        error_current = error_prop;
    end

    % Store chain and error
    chain(i,:) = params_current;
    errors(i) = error_current;
end
elapsed_time = toc;

%% === Select Samples After Burn-In ===
samples = chain(burn_in:end, :);  % select samples after burn-in
mean_estimates = mean(samples);

%% === Display Results ===
disp('Mean Parameter Estimates:');
disp(array2table(mean_estimates, 'VariableNames', param_names));

%% === Final Fit Plot ===
figure;
x_best = simulate_system(mean_estimates, x0, t, Vm);
plot(t, rad2deg(x(:,3)), 'b', 'DisplayName', 'True \alpha');
hold on;
plot(t, rad2deg(x_best(:,3)), 'r--', 'DisplayName', 'Estimated \alpha');
legend();
xlabel('Time (s)'); ylabel('Angle (deg)');
title('True vs Estimated Pendulum Response');
grid on;

%% === Error Plot (final) ===
figure;
semilogy(errors);
xlabel('Iteration'); ylabel('Error (log scale)');
title('Final Error Convergence');
grid on;

%% === Correlation Matrix & AIC ===
R = corrcoef(samples);
disp('Correlation Matrix:');
disp(array2table(R, 'VariableNames', param_names, 'RowNames', param_names));

best_error = compute_error(mean_estimates, x0, t, Vm, x);
AIC = 2 * param_dim + length(t) * log(best_error / length(t));
fprintf('AIC = %.2f\n', AIC);

%% === MCMC Parameter Estimation for Rotary Inverted Pendulum ===
% using the Metropolis-Hastings algorithm.
% This implementation follows the theoretical framework outlined in Chapter 9:
%   - Prior: Gaussian distribution on parameters
%   - Likelihood: Gaussian noise on measurements
%   - Posterior: log posterior = log likelihood + log prior
%   - Proposal: Gaussian random walk: cₜ = xₜ + ε,   ε ~ 𝒩(0, Σ)
%   - Proposal distribution: q(c|x) = 𝒩(x, Σ) = q(x|c) (symmetric)
%   - Acceptance: α(c,x) = min(1, π(c)/π(x))

clear; clc;

%% === Add Path and Load Simulated Data ===

addpath('D:\MATLAB\MCMC_Parameter_Estimation\Inverted_Pendulum_Model');
load('rip_sim_data.mat');  % loads: t, x

%% === Initial Setup ===

param_names = {'Br', 'Bp', 'kt', 'km', 'eta_m', 'eta_g'};
x0 = [0; 0; pi - 0.1; 0];  % initial state
Vm = @(t) chirp(t, 0.1, 100, 5, 'linear');  % input voltage

%% === Gaussian Priors (mean and std from similar motors) ===

% true_params = [0.0024, 0.0024, 0.007, 0.007, 0.69, 0.9];
mu_prior = [0.004, 0.004, 0.005, 0.005, 0.8, 0.8];     % prior mean µ
sigma_prior = [0.002, 0.002, 0.002, 0.002, 0.2, 0.2];  % prior std σ (uncertainty)

%% === MCMC Settings ===

num_iters = 10000;				% chain length N
param_dim = length(mu_prior);	% dimension d
burn_in = floor(num_iters / 2); % burn-in period

% === Proposal covariance Σ (diagonal, scaled by step size) ===
step_size = 0.01 * sigma_prior;

% === Parameter bounds ===
param_bounds = [1e-4, 1e-2;
                1e-4, 1e-2;
                1e-4, 1e-2;
                1e-4, 1e-2;
                0.5, 1.0;
                0.5, 1.0];

%% === Initialization ===
params_current = mu_prior;				% x₀ = prior mean
x_ref = x;								% reference data
chain = zeros(num_iters, param_dim);	% store all samples
errors = zeros(num_iters, 1);			% store errors

%% === Initial log posterior using compute_error ===

% MSE from model output
error_current = compute_error(params_current, x0, t, Vm, x_ref);	

noise_std = 0.05;% measurement noise stddev

% log-likelihood: ℓ(x) = -0.5 * (e/σ)²
log_like_current = -0.5 * (error_current / noise_std)^2;

% log-prior: log 𝒩(x|μ, σ²) = -0.5 * sum(((x - μ)/σ)²)
log_prior_current = -0.5 * sum(((params_current - mu_prior)./sigma_prior).^2);

% log-posterior: log π(x) = log ℓ(x) + log p(x)
log_post_current = log_like_current + log_prior_current;

%% === Initialize Real-Time Error Convergence Plot ===

figure('Name', 'MCMC Error Convergence', 'NumberTitle', 'off');
h_error = plot(1, errors(1), 'b');
xlabel('Iteration'); ylabel('Error');
title('Real-Time Error Convergence');
grid on;
hold on;

%% === Initialize Real-Time Parameter Trace Plots ===

figure('Name', 'Real-Time Parameter Trace', 'NumberTitle', 'off');
trace_plots = gobjects(param_dim, 1);
true_params = [0.0024, 0.0024, 0.007, 0.007, 0.69, 0.90];
for k = 1:param_dim
    subplot(param_dim, 1, k);
    trace_plots(k) = plot(1, params_current(k), 'b');
    % Plot true value (red horizontal line)
    yline(true_params(k), 'r--', 'LineWidth', 1.5, 'DisplayName', 'True Value');
	ylabel(param_names{k});
    grid on;
end
xlabel('Iteration');

fprintf('[%s] Starting MCMC Simulation...\n', datestr(now, 'HH:MM:SS'));

%% === MCMC Sampling Loop ===

tic;

for i = 1:num_iters
    fprintf('[%s] MCMC Iteration %d / %d\n', datestr(now, 'HH:MM:SS'), i, num_iters);

    % === Random Walk Proposal (Gaussian) ===
    % cₜ = xₜ + ε, ε ~ 𝒩(0, Σ)
    proposal = params_current + step_size .* randn(1, param_dim);

    % === Check bounds (hard constraint) === 
    if any(proposal < param_bounds(:,1)') || any(proposal > param_bounds(:,2)')
        chain(i,:) = params_current;
        errors(i) = error_current;
        continue;
    end

    % === Evaluate log-posterior of proposal === 

    % MSE from model output
	error_prop = compute_error(proposal, x0, t, Vm, x_ref);

    % log-likelihood: ℓ(x) = -0.5 * (e/σ)²
	log_like_prop = -0.5 * (error_prop / noise_std)^2;

    % log-prior: log 𝒩(x|μ, σ²) = -0.5 * sum(((x - μ)/σ)²)
	log_prior_prop = -0.5 * sum(((proposal - mu_prior)./sigma_prior).^2);

    % log-posterior: log π(x) = log ℓ(x) + log p(x)
	log_post_prop = log_like_prop + log_prior_prop;

    % === Compute Acceptance Probability α(c,x) ===
    % α(c,x) = min(1, π(c)/π(x))  ⇔  log α = log π(c) - log π(x)
    % Accept the proposal if log(U) < log α
    log_alpha = log_post_prop - log_post_current;  % log acceptance ratio

    if log(rand) < log_alpha
        % Proposal accepted: move to new sample c = proposal
        params_current = proposal;
        log_post_current = log_post_prop;
        error_current = error_prop;
    end

    % Store chain and error
    chain(i,:) = params_current;
    errors(i) = error_current;

	% Update real-time plots
	set(h_error, 'XData', 1:i, 'YData', errors(1:i));
	for k = 1:param_dim
    	set(trace_plots(k), 'XData', 1:i, 'YData', chain(1:i, k));
	end
	drawnow limitrate;

end
elapsed_time = toc;

%% === Post Processing - Select Samples After Burn-In ===

samples = chain(burn_in:end, :);  % select samples after burn-in
mean_estimates = mean(samples);

%% === Display Results ===

disp('Mean Parameter Estimates:');
disp(array2table(mean_estimates, 'VariableNames', param_names));

%% === Posterior Distributions ===

figure('Name', 'Posterior Distributions', 'NumberTitle', 'off');

for k = 1:param_dim
    subplot(param_dim, 1, k);
    histogram(samples(:,k), 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.4 0.6]);
    hold on;
    xline(mu_prior(k), '--r', 'LineWidth', 1.2, 'DisplayName', 'Prior Mean');
    xline(mean_estimates(k), '-k', 'LineWidth', 1.5, 'DisplayName', 'Posterior Mean');
    xlabel(param_names{k});
    ylabel('Density');
    legend();
    grid on;
end

sgtitle('Posterior Distributions of Parameters');

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
 
% MCMC Parameter Estimation for Rotary Inverted Pendulum
clear; clc;

% === Load Data ===
load('rip_sim_data.mat');  % contains: t, x

% === Initial Conditions ===
theta0  = 0; dtheta0 = 0;
alpha0  = pi - 0.1; dalpha0 = 0;
x0 = [theta0; dtheta0; alpha0; dalpha0];

% === Input Function ===
Vm = @(t) 1.5*sin(2*pi*0.5*t) + 1.0*sin(2*pi*1.5*t);

% === Parameters to Estimate ===
param_names = {'Br', 'Bp', 'kt', 'km'};

% === Prior Information ===
priors.mean = [0.0024, 0.0024, 0.007, 0.007];
priors.std  = [0.0005, 0.0005, 0.001, 0.001];
priors.lb   = [0.001, 0.001, 0.005, 0.005];
priors.ub   = [0.01, 0.01, 0.02, 0.02];

prior_type = 'exponential';  % Options: 'gaussian', 'lognormal', 'laplace', 'uniform', 'exponential'

% === MCMC Settings ===
num_iters = 5000;
burn_in = floor(num_iters / 2);
step_size = [1e-4, 1e-4, 1e-4, 1e-4];
params_curr = priors.mean;
chain = zeros(num_iters, 4);
error_curr = compute_error(params_curr, x0, t, Vm, x);

% === MCMC Loop ===
for i = 1:num_iters
    proposal = params_curr + step_size .* randn(1, 4);

    % Reject if out of bounds
    if any(proposal < priors.lb) || any(proposal > priors.ub)
        chain(i,:) = params_curr;
        continue;
    end

    error_prop = compute_error(proposal, x0, t, Vm, x);
    prior_curr = prior_density(params_curr, prior_type, priors);
    prior_prop = prior_density(proposal, prior_type, priors);

    post_curr = exp(-error_curr / 1e-3) * prior_curr;
    post_prop = exp(-error_prop / 1e-3) * prior_prop;

    alpha = min(1, post_prop / post_curr);

    if rand < alpha
        params_curr = proposal;
        error_curr = error_prop;
    end

    chain(i,:) = params_curr;
end

% === Analyze Samples ===
samples = chain(burn_in+1:end, :);
mean_est = mean(samples);
std_est  = std(samples);

disp('Estimated Parameters (mean ± std):');
for i = 1:4
    fprintf('%s: %.5f ± %.5f\n', param_names{i}, mean_est(i), std_est(i));
end

% === Plot Posterior Distributions ===
figure;
for i = 1:4
    subplot(2,2,i);
    histogram(samples(:,i), 40, 'Normalization', 'pdf');
    hold on;
    xline(priors.mean(i), '--r', 'Prior Mean');
    title(sprintf('Posterior of %s', param_names{i}));
    xlabel(param_names{i}); ylabel('Density');
    grid on;
end
sgtitle('Posterior Distributions');

% === Joint Posterior Plots ===
figure;
plotmatrix(samples);
sgtitle('Joint Posterior Distributions');

R = corrcoef(samples);
disp('Correlation Matrix:');
disp(array2table(R, 'VariableNames', param_names, 'RowNames', param_names));

best_error = compute_error(mean(samples), x0, t, Vm, x);
k = 4; n = length(t);
AIC = 2*k + n*log(best_error/n);
fprintf('AIC = %.2f\n', AIC);
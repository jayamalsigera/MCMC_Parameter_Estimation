clear; clc;

% === Load Simulated Data ===
load('rip_sim_data.mat');  % loads: t, x

% === Ground Truth Parameters ===
true_params = [0.0024, 0.0024, 0.007, 0.007];  % [Br, Bp, kt, km]

% === Initial Guess (Very Far) ===
initial_params = [0.001, 0.001, 0.001, 0.001];

% === Initial State ===
theta0  = 0; dtheta0 = 0;
alpha0  = pi - 0.1; dalpha0 = 0;
x0 = [theta0; dtheta0; alpha0; dalpha0];

% === Input Function ===
Vm = @(t) ...
    (t < 300) .* chirp(t, 0.1, 300, 5, 'linear') + ...
    (t >= 300 & t <= 600) .* chirp(t - 300, 0.1, 300, 5, 'logarithmic');

% === MCMC Settings ===
num_iters = 100000;
param_dim = length(initial_params);
step_size = [5e-4, 5e-4, 1e-3, 1e-3];  % Initial step size

param_bounds = [0.001, 0.01;   % Br
                0.001, 0.01;   % Bp
                0.005, 0.02;   % kt
                0.005, 0.02];  % km

params_current = initial_params;
chain = zeros(num_iters, param_dim);
errors = zeros(num_iters, 1);  % For convergence monitoring

% === Initial Error ===
error_prev = compute_error(params_current, x0, t, Vm, x);

% === MCMC Loop ===
for i = 1:num_iters
    % Adaptive step size adjustment every 2000 iterations (after burn-in)
    if mod(i, 2000) == 0 && i > 4000
        recent = chain(i-1999:i, :);
        recent_std = std(recent, 0, 1);
        step_size = 0.5 * recent_std + 1e-6;  % Add small offset to avoid zero
    end

    % Propose new sample
    proposal = params_current + step_size .* randn(1, param_dim);

    % Bounds check
    if any(proposal < param_bounds(:,1)') || any(proposal > param_bounds(:,2)')
        chain(i,:) = params_current;
        errors(i) = error_prev;
        continue;
    end

    % Simulate and compute error
    error_proposal = compute_error(proposal, x0, t, Vm, x);

    % Acceptance rule
    if error_proposal < error_prev || rand < exp((error_prev - error_proposal)/1e-3)
        params_current = proposal;
        error_prev = error_proposal;
    end

    % Store current state
    chain(i,:) = params_current;
    errors(i) = error_prev;
end

% === Plot Parameter Chains ===
figure;
for j = 1:param_dim
    subplot(param_dim,1,j);
    plot(chain(:,j));
    ylabel(sprintf('Param %d', j));
    title(sprintf('Chain for parameter %d', j));
end
xlabel('Iteration');

% === Mean Estimates After Burn-in ===
burn_in = floor(num_iters/2);
samples = chain(burn_in:end, :);
mean_estimates = mean(samples);

param_names = {'Br', 'Bp', 'kt', 'km'};
disp('Mean Parameter Estimates:');
disp(array2table(mean_estimates, 'VariableNames', param_names));

% === Posterior Distributions ===
figure;
for i = 1:4
    subplot(2,2,i);
    histogram(samples(:,i), 40, 'Normalization', 'pdf');
    hold on;
    xline(true_params(i), '--r', 'True');
    xlabel(param_names{i}); ylabel('Density');
    title(['Posterior of ', param_names{i}]);
    grid on;
end
sgtitle('Posterior Distributions');

% === Joint Posterior ===
figure;
plotmatrix(samples);
sgtitle('Joint Posterior Distributions');

% === Error Convergence Plot ===
figure;
semilogy(errors);
xlabel('Iteration'); ylabel('Error (log scale)');
title('Error Convergence Over Iterations');
grid on;

% === Correlation Matrix ===
R = corrcoef(samples);
disp('Correlation Matrix:');
disp(array2table(R, 'VariableNames', param_names, 'RowNames', param_names));

% === AIC Calculation ===
best_error = compute_error(mean_estimates, x0, t, Vm, x);
k = param_dim;
n = length(t);
AIC = 2*k + n*log(best_error/n);
fprintf('AIC = %.2f\n', AIC);
x_best = simulate_system(mean_estimates, x0, t, Vm);

figure;
plot(t, rad2deg(x(:,3)), 'r--', t, rad2deg(x_best(:,3)), 'b');
legend('True \alpha', 'Estimated \alpha');
xlabel('Time (s)'); ylabel('Angle (deg)');
title('True vs Estimated Pendulum Response');
grid on;
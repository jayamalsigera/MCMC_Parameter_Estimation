clear; clc;

% === Load Simulated Data ===
addpath('D:\MATLAB\MCMC_Parameter_Estimation\Inverted_Pendulum_Model');
load('rip_sim_data.mat');  % loads: t, x

% === Initial Setup ===
param_names = {'Br', 'Bp', 'kt', 'km'};
true_params = [0.0024, 0.0024, 0.007, 0.007];
initial_params = [0.001, 0.001, 0.001, 0.001];
x0 = [0; 0; pi - 0.1; 0];

% === Input Function (used during simulation) ===
Vm = @(t) chirp(t, 0.1, 100, 5, 'linear');

% === MCMC Parameters ===
num_iters = 10000;
param_dim = length(initial_params);
step_size = [1e-2, 1e-3, 1e-4, 1e-5];
% step_size = 1e-3;

param_bounds = [1e-3, 0.01;
                1e-3, 0.01;
                1e-3, 0.01;
                1e-3, 0.01];

params_current = initial_params;
chain = zeros(num_iters, param_dim);
errors = zeros(num_iters, 1);

% === Initial Error ===
error_prev = compute_error_mcmc(params_current, x0, t, Vm, x);

% === Initialize Real-Time Error Convergence Plot ===
figure('Name', 'MCMC Error Convergence', 'NumberTitle', 'off');
h_error = plot(1, error_prev, 'b');  % use plot() instead of semilogy
xlabel('Iteration'); ylabel('Error');
title('Real-Time Error Convergence');
grid on;
hold on;

% === Initialize Real-Time Parameter Trace Plots ===
figure('Name', 'Real-Time Parameter Trace', 'NumberTitle', 'off');
trace_plots = gobjects(param_dim, 1);  % Preallocate graphics objects

for k = 1:param_dim
    subplot(param_dim, 1, k);
    trace_plots(k) = plot(1, initial_params(k), 'b');
    ylabel(param_names{k});
    grid on;
end
xlabel('Iteration');
fprintf('[%s] Starting MCMC Simulation...\n', datestr(now, 'HH:MM:SS'));
tic;  % Start timing
% === MCMC Loop ===
for i = 1:num_iters
	if i == 1
    drawnow;  % force figure to render early
	end
	% === Update Real-Time Plots Every Iteration ===
	set(h_error, 'XData', 1:i, 'YData', errors(1:i));
	
	for k = 1:param_dim
    	set(trace_plots(k), 'XData', 1:i, 'YData', chain(1:i, k));
	end
	fprintf('[%s] MCMC Iteration %d / %d\n',datestr(now, 'HH:MM:SS'), i, num_iters);

	% === Adaptive step size based on 5 segments ===
	segment_length = floor(num_iters / 5);
	
	if mod(i, segment_length) == 0 && i >= segment_length
    	% Get recent segment
    	recent = chain(i - segment_length + 1 : i, :);
    	
    	% Compute standard deviation of parameter samples
    	recent_std = std(recent, 0, 1);
    	
    	% Update step size (scale down to encourage convergence)
    	step_size = 0.5 * recent_std + 1e-6;
	end

    % Propose new sample
    proposal = params_current + step_size .* randn(1, param_dim);

    % Reject if out of bounds
    if any(proposal < param_bounds(:,1)') || any(proposal > param_bounds(:,2)')
        chain(i,:) = params_current;
        errors(i) = error_prev;
        continue;
    end

    % Simulate proposal
    error_proposal = compute_error_mcmc(proposal, x0, t, Vm, x);

    % Acceptance rule
    if error_proposal < error_prev || rand < exp((error_prev - error_proposal)/1e-3)
        params_current = proposal;
        error_prev = error_proposal;
    end

    % Store
    chain(i,:) = params_current;
    errors(i) = error_prev;
	
	
	
drawnow limitrate;
	
end

elapsed_time = toc;  % End timing and store duration
fprintf('[%s] MCMC Simulation Completed!\n', datestr(now, 'HH:MM:SS'));
fprintf('MCMC Simulation Execution Time: %.2f seconds\n', elapsed_time);

% === Post-Processing ===
burn_in = floor(num_iters / 2);
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
    hold on; xline(true_params(i), '--r', 'True');
    title(['Posterior of ', param_names{i}]); xlabel(param_names{i}); ylabel('Density'); grid on;
end
sgtitle('Posterior Distributions');

% === Joint Posterior ===
figure;
plotmatrix(samples);
sgtitle('Joint Posterior Distributions');

% === Final Fit Plot ===
figure;
x_best = simulate_system_mcmc(mean_estimates, x0, t, Vm);
plot(t, rad2deg(x(:,3)), 'r--', t, rad2deg(x_best(:,3)), 'b');
legend('True \\alpha', 'Estimated \\alpha');
xlabel('Time (s)'); ylabel('Angle (deg)');
title('True vs Estimated Pendulum Response');
grid on;

% === Error Plot (final) ===
figure;
semilogy(errors);
xlabel('Iteration'); ylabel('Error (log scale)');
title('Final Error Convergence');
grid on;

% === Correlation Matrix & AIC ===
R = corrcoef(samples);
disp('Correlation Matrix:');
disp(array2table(R, 'VariableNames', param_names, 'RowNames', param_names));

best_error = compute_error_mcmc(mean_estimates, x0, t, Vm, x);
AIC = 2 * param_dim + length(t) * log(best_error / length(t));
fprintf('AIC = %.2f\\n', AIC);

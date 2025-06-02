function log_post = compute_log_posterior(params, x0, t, Vm, x_ref, mu_prior, sigma_prior)
    % === Log Likelihood ===
    x_sim = simulate_system_mcmc(params, x0, t, Vm);
    alpha_sim = x_sim(:,3);
    alpha_ref = x_ref(:,3);
    sigma_noise = 0.05;  % Assumed measurement noise
    log_likelihood = -0.5 * sum(((alpha_ref - alpha_sim) / sigma_noise).^2);

    % === Log Prior (Gaussian) ===
    log_prior = -0.5 * sum(((params - mu_prior) ./ sigma_prior).^2);

    % === Log Posterior ===
    log_post = log_likelihood + log_prior;
end

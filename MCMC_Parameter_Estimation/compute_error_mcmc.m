function err = compute_error_mcmc(params, x0, t, Vm, x_ref)
    x_sim = simulate_system_mcmc(params, x0, t, Vm);
    alpha_sim = x_sim(:,3);
    alpha_ref = x_ref(:,3);
    err = mean((alpha_sim - alpha_ref).^2);  % MSE
end

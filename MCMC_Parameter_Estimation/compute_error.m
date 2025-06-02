function err = compute_error(params, x0, t, Vm, x_ref)
    % COMPUTE_ERROR_MCMC Computes the weighted mean squared error (MSE) 
    % between simulated and reference system responses.
    %
    % Inputs:
    %   - params: parameter vector [Br, Bp, kt, km, eta_m, eta_g]
    %   - x0: initial state vector [theta; dtheta; alpha; dalpha]
    %   - t: time vector used for simulation
    %   - Vm: input voltage function handle, Vm(t)
    %   - x_ref: reference state trajectory (from measured or true data)
    %
    % Output:
    %   - err: scalar error value (weighted MSE)

    % === Simulate system with given parameters ===
    x_sim = simulate_system(params, x0, t, Vm);

    % === Define weights for each state variable in the error ===
    %   [theta, dtheta, alpha, dalpha] 
    %   alpha is more important, hence higher weight
    weights = [0, 0, 1, 1 ];

    % === Normalize weights so their sum is 1 ===
    weights = weights / sum(weights);

    % === Compute Weighted Mean Squared Error ===
    % Element-wise squared error * weights → sum across state variables → mean over time
    err = mean((x_sim - x_ref).^2 * weights');
end

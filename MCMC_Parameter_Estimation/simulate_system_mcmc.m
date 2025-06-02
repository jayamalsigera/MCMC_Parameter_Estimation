function x_sim = simulate_system(params, x0, t, Vm, iter)
    % SIMULATE_SYSTEM Integrates the RIP dynamics using given parameters.
    %
    % Inputs:
    %   - params: [Br, Bp, kt, km] (1x4 vector)
    %   - x0: initial state vector [theta; dtheta; alpha; dalpha]
    %   - t: time vector
    %   - Vm: input voltage function handle, Vm(t)
    %
    % Output:
    %   - x_sim: simulated state trajectory at specified time steps

    % === Unpack Parameters ===
    Br = params(1);
    Bp = params(2);
    kt = params(3);
    km = params(4);

    % === Integration with ode45 ===

    % Solve using ode45
    sol = ode45(@(t, x) rip_dynamics_mcmc(t, x, Vm(t), Br, Bp, kt, km), t, x0);

    % Evaluate solution at given time points
    x_sim = deval(sol, t)';  % Ensure output is Nx4 matrix
end

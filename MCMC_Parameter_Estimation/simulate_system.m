function x_sim = simulate_system(params, x0, t, Vm)
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
	Br     = params(1);
	Bp     = params(2);
	kt     = params(3);
	km     = params(4);
	eta_m  = params(5);
	eta_g  = params(6);


    % === Integration with ode45 ===

    % Solve using ode45
    sol = ode45(@(t, x) rip_dynamics(t, x, Vm(t), Br, Bp, kt, km, eta_m, eta_g), t, x0);

    % Evaluate solution at given time points
    x_sim = deval(sol, t)';  % Ensure output is Nx4 matrix
end

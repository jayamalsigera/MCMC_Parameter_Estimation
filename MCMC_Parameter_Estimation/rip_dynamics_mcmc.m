function dx = rip_dynamics_mcmc(t, x, Vm, Br, Bp, kt, km)
    % Extract states
    theta = x(1);
    dtheta = x(2);
    alpha = x(3);
    dalpha = x(4);

    % System constants (known, e.g., link lengths, masses)
    Jarm = 0.001;     % arm inertia
    Jpend = 0.001;    % pendulum inertia
    m = 0.2;          % pendulum mass
    l = 0.2;          % pendulum length
    g = 9.81;

    % Dynamics (example model — update to match your real model)
    ddtheta = (kt*Vm - Br*dtheta - m*l*g*sin(alpha)) / Jarm;
    ddalpha = (-Bp*dalpha + m*l*ddtheta*cos(alpha) - m*g*l*sin(alpha)) / Jpend;

    % State derivatives
    dx = zeros(4,1);
    dx(1) = dtheta;
    dx(2) = ddtheta;
    dx(3) = dalpha;
    dx(4) = ddalpha;
end

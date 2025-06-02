% === Initial Conditions ===
theta0  = 0;         % arm starts horizontal
dtheta0 = 0;         % arm at rest
alpha0  = pi - 0.1;  % pendulum starts near upright (pi) with small perturbation
dalpha0 = 0;         % pendulum at rest

x0 = [theta0; dtheta0; alpha0; dalpha0];

% === Input Voltage ===
% Vm = @(t) chirp(t, 0.1, 10, 3);  % From 0.1 Hz to 3 Hz in 10 s
Vm = @(t) ...
    (t < 300) .* chirp(t, 0.1, 300, 5, 'linear') + ...
    (t >= 300 & t <= 600) .* chirp(t - 300, 0.1, 300, 5, 'logarithmic');

% === ODE Solver Options ===
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

% === Simulation ===
tspan = [0 600];
[t, x] = ode45(@(t, x) rip_dynamics(t, x, Vm(t)), tspan, x0, opts);

% === Save as .mat file ===
save('rip_sim_data.mat', 't', 'x');

% === Plot in Degrees ===
theta_deg = rad2deg(x(:,1));
alpha_deg = rad2deg(x(:,3));

% === Wrap Angles to [-180, +180] ===
theta_deg_wrapped = mod(theta_deg + 180, 360) - 180;
alpha_deg_wrapped = mod(alpha_deg + 180, 360) - 180;

% === Plot ===
figure;
plot(t, theta_deg_wrapped, 'b', t, alpha_deg_wrapped, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle (°)');
legend('\theta (arm)', '\alpha (pendulum)');
title('Rotary Inverted Pendulum Response (Degrees, Wrapped)');
grid on;
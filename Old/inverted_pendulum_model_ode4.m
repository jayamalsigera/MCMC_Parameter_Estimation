% === Initial Conditions ===
theta0  = 0;
dtheta0 = 0;
alpha0  = pi - 0.1;
dalpha0 = 0;
x0 = [theta0; dtheta0; alpha0; dalpha0];

% === Input Voltage ===
Vm = @(t) ...
    (t < 300) .* chirp(t, 0.1, 300, 5, 'linear') + ...
    (t >= 300 & t <= 600) .* chirp(t - 300, 0.1, 300, 5, 'logarithmic');

% === Time Vector for Fixed-Step ===
dt = 0.001;
tspan = 0:dt:600;
N = length(tspan);
x = zeros(N, 4);
x(1, :) = x0';

% === RK4 Fixed-Step Integration ===
for i = 1:N-1
    t = tspan(i);
    u = Vm(t);
    
    k1 = rip_dynamics(t, x(i,:)', u);
    k2 = rip_dynamics(t + dt/2, x(i,:)' + dt/2 * k1, u);
    k3 = rip_dynamics(t + dt/2, x(i,:)' + dt/2 * k2, u);
    k4 = rip_dynamics(t + dt, x(i,:)' + dt * k3, u);
    
    x(i+1, :) = x(i,:) + (dt/6)*(k1' + 2*k2' + 2*k3' + k4');
end

% === Save as .mat file ===
save('rip_sim_data.mat', 'tspan', 'x');

% === Plot in Degrees ===
theta_deg = rad2deg(x(:,1));
alpha_deg = rad2deg(x(:,3));
theta_wrapped = mod(theta_deg + 180, 360) - 180;
alpha_wrapped = mod(alpha_deg + 180, 360) - 180;

figure;
plot(tspan, theta_wrapped, 'b', tspan, alpha_wrapped, 'r--', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Angle (°)');
legend('\theta (arm)', '\alpha (pendulum)');
title('RIP Response (Fixed Step 0.001s, Degrees, Wrapped)');
grid on;

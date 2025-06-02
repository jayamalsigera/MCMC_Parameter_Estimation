function dx = rip_dynamics(t, x, Vm)
% Rotary Inverted Pendulum Nonlinear Dynamics
% x = [theta; dtheta; alpha; dalpha]

% === System Parameters ===

% Pendulum
mp  = 0.127;
Lp  = 0.337;
lp  = 0.156;
Jp  = 0.0012 + mp*lp^2;
Bp  = 0.0024;

% Rotary arm
mr  = 0.257;
Lr  = 0.216;
Jr  = 0.0020;
Br  = 0.0024;

% Motor
Kg    = 70;
km    = 0.007;
kt    = 0.007;
Rm    = 2.6;
eta_m = 0.69;
eta_g = 0.90;

g = 9.81;

% === States ===
theta  = x(1);
dtheta = x(2);
alpha  = x(3);
dalpha = x(4);

% === Motor torque ===
tau = (eta_g * Kg * eta_m * kt * (Vm - Kg * km * dtheta)) / Rm;

% === Dynamics ===
c = cos(alpha);
s = sin(alpha);

M11 = Jr + mp*Lr^2 + mp*lp^2 - mp*lp^2*c^2;
M12 = -mp * Lr * lp * c;
M21 = M12;
M22 = Jp + mp * lp^2;
M = [M11, M12; M21, M22];

RHS1 = tau - Br*dtheta + (mp*lp^2*s*c + mp*Lr*lp*s)*dalpha*dtheta;
RHS2 = -Bp*dalpha - mp*lp^2*c*s*dtheta^2 - mp*lp*g*s;
RHS = [RHS1; RHS2];

acc = M \ RHS;

dx = zeros(4,1);
dx(1) = dtheta;
dx(2) = acc(1);
dx(3) = dalpha;
dx(4) = acc(2);
end

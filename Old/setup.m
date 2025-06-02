% === Parameters ===

% Pendulum (vertical link)
mp  = 0.127;             % kg
Lp  = 0.337;             % m
lp  = 0.156;             % m (CoM)
Jp  = 0.0012 + mp*lp^2;  % kg*m^2 (about pivot using parallel axis)
Bp  = 0.0024;            % N*m*s/rad

% Rotary arm (horizontal link)
mr  = 0.257;             % kg
Lr  = 0.216;             % m
Jr  = 0.0020;            % kg*m^2 (about pivot)
Br  = 0.0024;            % N*m*s/rad

% Motor and gear parameters
Kg     = 70;             % gear ratio
km     = 0.007;          % V*s/rad (back EMF)
kt     = 0.007;          % N*m/A
Rm     = 2.6;            % Ohms
eta_m  = 0.69;           % motor efficiency
eta_g  = 0.90;           % gear efficiency

% Gravity
g = 9.81;                % m/s^2



% % Eq. 1
% left1 = (mp*Lr^2 + 0.25*mp*Lp^2 - mp*Lp^2*cos(alpha)^2 + Jr) * ddtheta;
% left2 = -0.5 * mp * Lp * Lr * cos(alpha) * ddalpha;
% left3 = (0.5 * mp * Lp^2 * sin(alpha) * cos(alpha) + 0.5 * mp * Lp * Lr * sin(alpha)) * dalpha * dtheta;
% right1 = tau - Br * dtheta;
% 
% % Eq. 2
% left4 = -0.5 * mp * Lp * Lr * cos(alpha) * ddtheta;
% left5 = (Jp + 0.25 * mp * Lp^2) * ddalpha;
% left6 = -0.25 * mp * Lp^2 * cos(alpha) * sin(alpha) * dtheta^2;
% right2 = -0.5 * mp * Lp * g * sin(alpha) - Bp * dalpha;
% 
% % Eq. 3
% tau = eta_g * Kg * eta_m * kt * (Vm - Kg * km * dtheta) / Rm;
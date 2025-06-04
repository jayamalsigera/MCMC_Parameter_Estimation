function dx = rip_dynamics(t, x, Vm, Br, km, eta_m, eta_g)
    % ROTARY INVERTED PENDULUM NONLINEAR DYNAMICS (for MCMC Estimation)
    % 
    % Inputs:
    %   t       - current simulation time (not used here, but required by ODE)
    %   x       - state vector: [theta; dtheta; alpha; dalpha]
    %   Vm      - input motor voltage (scalar)
    %   Br      - arm damping coefficient (to estimate)
    %   Bp      - pendulum damping coefficient (to estimate)
    %   kt      - motor torque constant (to estimate)
    %   km      - back-emf constant (to estimate)
    %   eta_m   - motor efficiency (to estimate)
    %   eta_g   - gear efficiency (to estimate)
    %
    % Output:
    %   dx      - time derivative of the state vector

    % === Physical Constants ===
    
    % Pendulum properties
    mp  = 0.127;      % Pendulum mass (kg)
    Lp  = 0.337;      % Full pendulum length (m)
    lp  = 0.156;      % Distance to pendulum center of mass (m)
    Jp  = 0.0012 + mp*lp^2;  % Pendulum moment of inertia (kg.m^2)
    
    % Rotary arm properties
    mr  = 0.257;      % Arm mass (kg)
    Lr  = 0.216;      % Arm length (m)
    Jr  = 0.0020;     % Arm moment of inertia (kg.m^2)
    
    % Motor and gear properties
    Kg    = 70;       % Gear ratio
    Rm    = 2.6;      % Motor resistance (Ohms)
    g     = 9.81;     % Gravity (m/s^2)


	% === ESTIMATION PARAMETERS ===
	% Bp
	% 
	% km
	% kt
	% eta_m
	% eta_g
	Bp = Br;
	kt = km;

    % === Extract States ===
    theta  = x(1);    % Arm angle (rad)
    dtheta = x(2);    % Arm angular velocity (rad/s)
    alpha  = x(3);    % Pendulum angle (rad)
    dalpha = x(4);    % Pendulum angular velocity (rad/s)

    % === Motor Torque Equation ===
    % tau = ηg * Kg * ηm * kt * (Vm - Kg * km * dtheta) / Rm
    tau = (eta_g * Kg * eta_m * kt * (Vm - Kg * km * dtheta)) / Rm;

    % === Precompute Trig Functions ===
    c = cos(alpha);
    s = sin(alpha);

    % === Mass Matrix ===
    % M * [ddtheta; ddalpha] = RHS
    M11 = Jr + mp*Lr^2 + mp*lp^2 - mp*lp^2*c^2;
    M12 = -mp * Lr * lp * c;
    M21 = M12;
    M22 = Jp + mp * lp^2;
    M = [M11, M12; M21, M22];

    % === Right-Hand Side (Forces and Torques) ===
    RHS1 = tau ...
         - Br*dtheta ...
         + (mp*lp^2*s*c + mp*Lr*lp*s) * dalpha * dtheta;

    RHS2 = -Bp*dalpha ...
         - mp*lp^2*c*s*dtheta^2 ...
         - mp*lp*g*s;

    RHS = [RHS1; RHS2];

    % === Solve for Angular Accelerations ===
    acc = M \ RHS;

    % === Return State Derivative ===
    dx = zeros(4,1);
    dx(1) = dtheta;     % d(theta)/dt = angular velocity
    dx(2) = acc(1);     % d(dtheta)/dt = arm angular acceleration
    dx(3) = dalpha;     % d(alpha)/dt = pendulum angular velocity
    dx(4) = acc(2);     % d(dalpha)/dt = pendulum angular acceleration
end

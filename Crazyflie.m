% =========================================================================
% Crazyflie 2.1 Nonlinear Model Simulation
% =========================================================================
clear; clc; close all;

% Crazyflie 2.1 Parameters 
m = 0.027; % Mass in kg
g = 9.81;  % Gravity (m/s^2)
l = 0.05;  % Arm length in m
% Diagonal Inertia Matrix J = diag(Jx, Jy, Jz) in kg*m^2
Jx = 3.3e-5; Jy = 3.6e-5; Jz = 5.9e-5; 
J = diag([Jx, Jy, Jz]); 


%% Simulation
dt=0.01;
t = 0:dt:5; % Simulation interval
X = zeros(12, length(t));

% Initial state vector: [p; v; lambda; omega]
% p=[14;15;16], v=[0;0;0], lambda=[0;0;0], omega=[0;0;0]  
X(:, 1) = zeros(12, 1);
X(1:3)=[14 15 16];

% Inputs: Small constant thrust, equal to all rotors
% If T < m*g, it stays on the ground. To hover apply thrust higher than
% drone wheight
% T/rotor = c ; Total T = 4*c.
T_total = m * g; 
np = [0; 0; 0]; % Equal thrust means zero moments

for i = 1:(length(t) - 1)
    current_t = t(i);
    current_X = X(:, i);

    dX = drone_dynamics(current_t, current_X, m, g, J, T_total, np);

    %    Next State = Current State + (Rate of Change * Time Step)
    X(:, i+1) = current_X + (dX * dt);
end


% 5. Plot Results (Z-Position and Z-Velocity)
figure;
subplot(2,1,1);
plot(t, X(3,:), 'LineWidth', 2);
title('Crazyflie Z-Position (Altitude)');
xlabel('Time (s)'); ylabel('z (m)');
grid on;

subplot(2,1,2);
plot(t, X(6,:), 'LineWidth', 2);
title('Crazyflie Z-Velocity');
xlabel('Time (s)'); ylabel('w (m/s)');
grid on;

% O ~ é só pq se passa o tempo para dntro da função mas não é utilizado
function dX = drone_dynamics(~, X, m, g, J, T, np)
    % Group state variables
    p = X(1:3);
    v = X(4:6);
    lambda = X(7:9);
    phi = lambda(1); theta = lambda(2); psi = lambda(3);
    omega = X(10:12);
    
    % Rotation Matrix (Z-Y-X Euler to Body->Inertial)
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    R = Rz * Ry * Rx; % Body to Inertial
    
    % Euler Rate Matrix
    Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
     
    % Gravity Force in Inertial Frame (ENU convention => Z up)
    fg_inertial = [0; 0; -m*g];
    
    % Drag matrix
    fa = [0; 0; 0];
    
    % 1. Position Kinematics (Body velocity to Inertial velocity)
    dp = R * v;
    % 2. Velocity Dynamics (Newton's 2nd Law in Body Frame)
    % R' * fg_inertial rotates gravity into the body frame
    dv = cross(-omega, v) + (1/m) * (R' * fg_inertial + fa + [0; 0; T]);
    % 3. Attitude Kinematics
    dlambda = Q * omega;
    % 4. Angular Velocity Dynamics (Euler's Equations)
    domega = J \ (cross(-omega, J * omega) + np);
    
    % Derivatives
    dX = [dp; dv; dlambda; domega];
end
% =========================================================================
% Crazyflie Drone motion control and ICUAS competition
% =========================================================================
clear; clc; close all;

%% 0. Global Symbolic Parameters and Operating Point
% Symbolic state and input
syms px py pz vx vy vz real 
syms uax uay uaz real       
X_sym = [px; py; pz; vx; vy; vz];
u_a_sym = [uax; uay; uaz];

% Symbolic parameters and operation point variables
syms m g l dx dy dz T phi theta psi real
param_syms = [m, g, dx, dy, dz, T, phi, theta, psi];

% Numeric values for the physical drone
m_num = 0.027; % Mass in kg
g_num = 9.81; % Gravity in m/s^2
l_num = 0.05; % Arm length in m
dx_num = 9.18e-7;      
dy_num = 9.18e-7;      
dz_num = 10.31e-7;

% Operation Point (e.g., Hover)
T_op = m_num * g_num; % Hover thrust compensates gravity
phi_op = 0; % 0 radians roll
theta_op = 0; % 0 radians pitch
psi_op = 0; % 0 radians yaw

% Group numeric values corresponding to param_syms
param_values = [m_num, g_num, dx_num, dy_num, dz_num, T_op, phi_op, theta_op, psi_op];

%% 1.1 Simplified Nonlinear model 
% (Function drone_nonlinear_dynamics is defined at the bottom)

%% 1.2 Small roll and pitch angles approximation, zero yaw
D_sym = diag([dx, dy, dz]); 

% Symbolic u_lambda
u_lambda_sym = [T, phi, theta, psi];

% Symbolic rotation and equivalent acceleration
R_small_sym = small_angles_0_yaw_R(u_lambda_sym);
u_a_eq_sym  = u_lamdaTOu_a(u_lambda_sym, m, g);

%% 1.3 Define the linear model using symbolic Jacobians
% Symbolic dynamics
dX_sym = drone_small_angle_dynamics(X_sym, u_a_sym, m, D_sym, R_small_sym);

% Symbolic Jacobians
A_sym = jacobian(dX_sym, X_sym);
B_sym = jacobian(dX_sym, u_a_sym);

% Substitute the global numeric values into the symbolic Jacobians
A = double(subs(A_sym, param_syms, param_values));
B = double(subs(B_sym, param_syms, param_values));

%% 1.4 Design an LQR controller
qp = 50;
qv = 1;
r  = 10;

Q = [qp*eye(3),  zeros(3,3);
     zeros(3,3), qv*eye(3)];
R_lqr = r * eye(3);

[K, P, poles] = lqr(A, B, Q, R_lqr);

disp('--- LQR Gain Matrix K ---');
disp(K);

% Uncomment to check closed-loop stability
% eig_drone_LQR = eig(A - B*K);
% disp('Closed-loop poles:'); disp(eig_drone_LQR);


%% ========================================================================
% Functions
% =========================================================================

%--------------------------------------------------------------------------
% For 1.1: Full Nonlinear Dynamics
function dX = drone_nonlinear_dynamics(X, u_lambda, m, g, D)
    p = X(1:3);
    v = X(4:6);
    T = u_lambda(1);
    phi = u_lambda(2); 
    theta = u_lambda(3); 
    psi = u_lambda(4);
    
    % Z-Y-X Euler Rotation
    R = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
         -sin(theta)        , sin(phi)*cos(theta)                             , cos(phi)*cos(theta)                              ];
         
    dp = v; 
    dv = ([0;0;-m*g] - R*D*R'*v + R*[0;0;T]) / m;
    
    dX = [dp; dv];
end

%--------------------------------------------------------------------------
% For 1.2: Input Transformations
function u_a = u_lamdaTOu_a(u_lambda, m, g)
    T = u_lambda(1);
    phi = u_lambda(2);
    theta = u_lambda(3);
    
    u_a = [T*theta/m; -T*phi/m; T/m - g];
end

function R_small = small_angles_0_yaw_R(u_lambda)
    phi = u_lambda(2);
    theta = u_lambda(3);
    % Confirmar q esta matriz está corrreta
    R_small = [1, phi*theta, theta;
               0, 1, -phi;
              -theta, phi, 1];
end

%--------------------------------------------------------------------------
% For 1.3: Linearized Dynamics
function dX = drone_small_angle_dynamics(X, u_a, m, D, R)
    p = X(1:3);
    v = X(4:6);
    
    dp = v; 
   
    dv = u_a - (R * D * R' * v) / m;
    
    dX = [dp; dv];
end
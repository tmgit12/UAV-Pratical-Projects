% =========================================================================
% Crazyflie Drone Modeling and Identification
% =========================================================================
clear; clc; close all;

% Crazyflie 2.1 Parameters 
m = 0.027; % Mass in kg
g = 9.81; % Gravity in m/s^2
l = 0.05; % Arm length in m
Jx = 3.3e-5; Jy = 3.6e-5; Jz = 5.9e-5; 
J = diag([Jx, Jy, Jz]); % Inertia Matrix in kg*m^2
dx = 0.01; dy = 0.01; dz = 0.01;
D = diag([dx, dy, dz]); % Drag Matrix in kg/s

%% 1.5 Nonlinear model simulation
dt=0.01;
t = 0:dt:5; % Simulation interval
X = zeros(12, length(t));
U = zeros(4,length(t));  % U = [T; n_p] (4 inputs: thrust + 3 moments)

% Initial state vector: [p; v; lambda; omega]
% p=[14;15;16], v=[0;0;0], lambda=[0;0;0], omega=[0;0;0]  
X(1:3, 1) = [14; 15; 16];

% Inputs: Small constant thrust, equal to all rotors
% T_Total = T1+T2+T3+T4
T_total = 4*m * g; 
n_p = [0; 0; 0]; % Equal thrust means zero moments
U = [T_total; n_p] * ones(1, length(t));



for i = 1:(length(t) - 1)
    current_t = t(i);
    current_X = X(:, i);
    current_U = U(:, i);

    dX = drone_dynamics(current_X, m, g, J, current_U, D);

    % Next State = Current State + (Rate of Change * Time Step)
    X(:, i+1) = current_X + (dX * dt);
end

%% Plotting
% 3D Trajectory
figure('Name', '3D Trajectory');
h1 = plot3(X(1,:), X(2,:), X(3,:), 'b-', 'LineWidth', 1.5);
grid on; hold on;
h2 = plot3(X(1,1), X(2,1), X(3,1), 'go', 'MarkerFaceColor', 'g'); 
h3 = plot3(X(1,end), X(2,end), X(3,end), 'ro', 'MarkerFaceColor', 'r'); 
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
legend([h1, h2, h3], {'Trajectory', 'Initial', 'Final'});
title('Flight Path');

% Euler Angles
figure('Name', 'Attitude');
plot(t, X(7:9,:) * (180/pi));
xlabel('t [s]');
ylabel('Attitude [deg]');
grid on;
legend('\phi (roll)', '\theta (pitch)', '\psi (yaw)');

% Velocity
figure('Name', 'Velocity');
plot(t, X(4:6,:));
xlabel('t [s]');
ylabel('Velocity [m/s]');
grid on;
legend('u', 'v', 'w');

%% 1.6 Determining x_e and u_e for OP2
% Defining vx_e = 0; vy_e = 0; falls under the conditions of OP1. OP1 is a
% particularity of OP2.

syms m g Jx Jy Jz dx dy dz real % Crazyflie parameters

% Desired horizontal equilibrium trajectory:
Vx_e = 0.5; Vy_e = 0; psi_e = pi/4;
% Here Vx_e and Vy_e are the velocities in the inertial frame represented
% in the report as \dot{x}_0 and \dot{y}_0. So Vx_e, Vy_e and psi_e
% represent here the desired horizontal equilibrium trajectory labeled in
% the report as \dot{x}_0, \dot{y}_0 and \psi_0 respectively.

% States and inputs as functions of the desired equilibrium parameters:
[u_e, v_e, w_e, theta_e, phi_e, T_e] = State_input_equilibrium_functions(Vx_e, Vy_e, psi_e, m, g, dx, dy, dz);

%% 1.7 Linearized models

syms px py pz vx vy vz phi theta psi wx wy wz real % States
% We will use this notation instead of X=[x,y,z, u,v,w, phi,theta,psi, p,q,
% r] to avoid conflict with other variables 
syms T nx ny nz real % Inputs

% Grouping:
X = [px; py; pz; vx; vy; vz; phi; theta; psi; wx; wy; wz];
U = [T; nx; ny; nz];
J = diag([Jx, Jy, Jz]);
D = diag([dx, dy, dz]);

% dX=f(x,u) symbolic dynamics
f = drone_dynamics(X, m, g, J, U, D);

% Symbolic Jacobians
A_sym = jacobian(f, X);
B_sym = jacobian(f, U);

% Physical parameter value substitution
param_syms   = [m, g, Jx, Jy, Jz, dx, dy, dz];
param_values = [0.027, 9.81, 3.3e-5, 3.6e-5, 5.9e-5, 0.01, 0.01, 0.01];

%--------------------------------------------------------------------------
% OP1 (Hover)
Vx_op1 = 0; Vy_op1 = 0; psi_op1 = 0; % Stationary hover
[u1, v1, w1, th1, ph1, T1] = State_input_equilibrium_functions(Vx_op1, Vy_op1, psi_op1, ...
                                    param_values(1), param_values(2), param_values(6), ...
                                    param_values(7), param_values(8));

% X_e_OP1=[p; v; lambda; omega]
% Position is arbitrary, we'll use [14;15;16] Angular velocities are 0.
xe_op1 = [14; 15; 16; u1; v1; w1; ph1; th1; psi_op1; 0; 0; 0];
ue_op1 = [T1; 0; 0; 0];

% Substitute parameters and equilibrium values to get numeric A and B
A_OP1 = double(subs(A_sym, [X; U; param_syms'], [xe_op1; ue_op1; param_values']));
B_OP1 = double(subs(B_sym, [X; U; param_syms'], [xe_op1; ue_op1; param_values']));

%--------------------------------------------------------------------------
% OP2 (Horizontal Flight)
Vx_op2 = 0.5; Vy_op2 = 0; psi_op2 = 10*pi/180; % 0.5 m/s forward flight at 45 deg yaw
[u2, v2, w2, th2, ph2, T2] = State_input_equilibrium_functions(Vx_op2, Vy_op2, psi_op2, ...
                                    param_values(1), param_values(2), param_values(6), ...
                                    param_values(7), param_values(8));

% X_e_OP1=[p; v; lambda; omega]
% Position is arbitrary, we'll use [14;15;16] Angular velocities are 0.
xe_op2 = [14; 15; 16; u2; v2; w2; ph2; th2; psi_op2; 0; 0; 0];
ue_op2 = [T2; 0; 0; 0];

A_OP2 = double(subs(A_sym, [X; U; param_syms'], [xe_op2; ue_op2; param_values']));
B_OP2 = double(subs(B_sym, [X; U; param_syms'], [xe_op2; ue_op2; param_values']));

%% 1.8 Eigenvalues
eig_OP1 = eig(A_OP1);
eig_OP2 = eig(A_OP2);

figure('Name', 'Eigenvalues Plot');
plot(real(eig_OP1), imag(eig_OP1), 'bx', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(real(eig_OP2), imag(eig_OP2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(0, 'k--'); yline(0, 'k--');
grid on;
xlabel('Real'); ylabel('Imaginary');
title('Eigenvalues: OP1 (Hover) vs OP2 (Horizontal)');
legend('OP1 (Hover)', 'OP2 (Horizontal)');

%% 1.9 Analyze controlability, observability and stability

%--------------------------------------------------------------------------
% OP1

n_OP1 = length(A_OP1);
C_OP1 = eye(12);

% Controlability
Mc_OP1 = ctrb( A_OP1, B_OP1);  % Mc = [ B | A*B | ... | A^(n-1)*B ]
rank(Mc_OP1);
if ( rank(Mc_OP1) == n_OP1 )
    disp( 'OP1 linear model is fully controllable' )
else
    disp( 'OP1 linear model is not fully controllable')
end

% Observability
Mo_OP1 = obsv(A_OP1 , C_OP1);  % Mo = [ C; C*A; ...; C*A^(n-1) ]

rank(Mo_OP1);
if ( rank(Mo_OP1) == n_OP1 )
    disp( 'OP1 linear model is fully observable' )
else
    disp( 'OP1 linear model is not fully observable')
end

% Stability
J_OP1 = jordan(A_OP1);

%--------------------------------------------------------------------------
% OP2

n_OP2 = length(A_OP2);
C_OP2 = eye(12);

% Controlability
Mc_OP2 = ctrb( A_OP2, B_OP2);  % Mc = [ B | A*B | ... | A^(n-1)*B ]
rank(Mc_OP2);
if ( rank(Mc_OP2) == n_OP2 )
    disp( 'OP2 linear model is fully controllable' )
else
    disp( 'OP2 linear model is not fully controllable')
end

% Observability
Mo_OP2 = obsv(A_OP2 , C_OP2);  % Mo = [ C; C*A; ...; C*A^(n-1) ]

rank(Mo_OP2);
if ( rank(Mo_OP2) == n_OP2 )
    disp( 'OP2 linear model is fully observable' )
else
    disp( 'OP2 linear model is not fully observable')
end

% Stability
J_OP2 = jordan(A_OP2);

%% 1.10 Deduce transfer functions
%--------------------------------------------------------------------------
% OP1

sys_OP1 = ss(A_OP1, B_OP1, eye(12), zeros(12,4));
sys_tf_OP1 = tf(sys_OP1);

% Inner Loops (Attitude) - Use zpk to see the exact poles and zeros
G_nx_phi_sym_OP1 = zpk(minreal(sys_tf_OP1(7, 2))); % nx/phi
G_ny_theta_OP1 = zpk(minreal(sys_tf_OP1(8, 3))); % ny/theta
G_nz_psi_OP1 = zpk(minreal(sys_tf_OP1(9, 4))); % nz/psi

% Outer Loops (Position) - minreal
G_px_theta_OP1 = zpk(minreal(sys_tf_OP1(1, 3) / sys_tf_OP1(8, 3))); % px/theta
G_py_phi_OP1 = zpk(minreal(sys_tf_OP1(2, 2) / sys_tf_OP1(7, 2))); % py/phi
G_pz_T_OP1 = zpk(minreal(sys_tf_OP1(3, 1))); % pz/T

%--------------------------------------------------------------------------
% OP2
sys_OP2 = ss(A_OP2, B_OP2, eye(12), zeros(12,4));
sys_tf_OP2 = tf(sys_OP2);

% Inner Loops (Attitude) - Use zpk to see the exact poles and zeros
G_nx_phi_sym_OP2 = zpk(minreal(sys_tf_OP2(7, 2))); % nx/phi
G_ny_theta_OP2 = zpk(minreal(sys_tf_OP2(8, 3))); % ny/theta
G_nz_psi_OP2 = zpk(minreal(sys_tf_OP2(9, 4))); % nz/psi

% Outer Loops (Position) - minreal
G_px_theta_OP2 = zpk(minreal(sys_tf_OP2(1, 3) / sys_tf_OP2(8, 3))); % px/theta
G_py_phi_OP2 = zpk(minreal(sys_tf_OP2(2, 2) / sys_tf_OP2(7, 2))); % py/phi
G_pz_T_OP2 = zpk(minreal(sys_tf_OP2(3, 1))); % pz/T

%% 1.11 Closed loop stability of G_{T,pz}(s) and G_{ϕ,py}(s)

%--------------------------------------------------------------------------
% OP1

figure;
rlocus(G_pz_T_OP1);
figure;
rlocus(G_py_phi_OP1);

%--------------------------------------------------------------------------
% OP2
figure;
rlocus(G_pz_T_sym_OP2);
figure;
rlocus(G_py_phi_sym_OP2);

%% 2. Identification
clear all

% read file
csvfilename = '2024-04-04_log07.csv';

array = dlmread(csvfilename,',',1,0);

% get data from table (octave)
time = array(:,1)'*1e-3;
%dt = mean(t(2:end)-t(1:end-1)); % estimate of sample time
dt = mean(diff(time));
pos = array(:,2:4)';
vel = array(:,5:7)';
lbd = array(:,8:10)'*pi/180;
om = array(:,11:13)'*pi/180;
pos_ref = array(:,14:16)';
yaw_ref = array(:,17)';
motors = array(:,18:21)';

% convert date to print format
t = time - time(1);
x = [pos;vel;lbd;om];
x_ref = [pos_ref;0*vel;lbd*0;om*0];
x_ref(9,:) = yaw_ref;
uint16_max = 2^16;
u = motors/uint16_max;

% plot data
%initPlots;
%vehicle3d_ref_show_data(t,x,u,x_ref);

% prepare data for ID
% dt = diff(t);
% dt = [dt,dt(end)];
% sample_time_stats = [mean(dt),min(dt),max(dt)],
% Ts = sample_time_stats(1);
%u_id = lbd(1,:)';
%y_id = pos(2,:)';

%% 2.1 Motion in x,y and z

% Looking at the position plots the following intervals of time were
% considered for movement in each axis:
t_x_start = 74;  t_x_end = 101; % Replace with actual times from plot
t_y_start = 104; t_y_end = 131; % Replace with actual times from plot
t_z_start = 44; t_z_end = 71; % Replace with actual times from plot

% indices dos intervalos de tempo
idx_X = (t >= t_x_start) & (t <= t_x_end);
idx_Y = (t >= t_y_start) & (t <= t_y_end);
idx_Z = (t >= t_z_start) & (t <= t_z_end);

px = pos(1,:);
py = pos(2,:);
pz = pos(3,:);

phi = lbd(1,:);
theta = lbd(2,:);
psi = lbd(3,:);

% Datasets para movimento no x, y e z
Dataset_X.t = t(idx_X) - t(find(idx_X, 1, 'first'));
Dataset_X.px = px(idx_X);
Dataset_X.theta = theta(idx_X);

Dataset_Y.t = t(idx_Y) - t(find(idx_Y, 1, 'first'));
Dataset_Y.py = py(idx_Y);
Dataset_Y.phi = phi(idx_Y);

Dataset_Z.t = t(idx_Z) - t(find(idx_Z, 1, 'first'));
Dataset_Z.pz = pz(idx_Z);
Dataset_Z.thrust = u(:, idx_Z);

%% 2.3 Identify modes

data_x = iddata(Dataset_X.px', Dataset_X.theta', dt);
data_y = iddata(Dataset_Y.py', Dataset_Y.phi', dt);
data_z = iddata(Dataset_Z.pz', Dataset_Z.thrust', dt);

% Change the last two values corresponding to the expected amount of zeros
% and poles, respectively, for each transfer function
sys_px_theta = tfest(detrend(data_x), 2, 0);
sys_py_phi = tfest(detrend(data_y), 2, 0);
sys_pz_T = tfest(detrend(data_z), 2, 0);

%% 2.4 Model comparison

figure;
compare(data_x, sys_px_theta);

figure;
compare(data_y, sys_py_phi);

figure;
compare(data_z, sys_pz_T);

%% Functions

function S = skew(x)
    S = [ 0    -x(3)  x(2);
          x(3)  0     -x(1);
         -x(2)  x(1)   0 ];
end

function dX = drone_dynamics(X, m, g, J, u, D)
    % Group state variables
    p = X(1:3);
    v = X(4:6);
    lambda = X(7:9);
    phi = lambda(1); theta = lambda(2); psi = lambda(3);
    omega = X(10:12);

    T=u(1);
    np=u(2:4);

    % Skew-symmetric matrix for omega
    S_omega = skew(omega);
    
    % Rotation Matrix (Z-Y-X Euler to Body->Inertial)
    R = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
         -sin(theta)        , sin(phi)*cos(theta)                             , cos(phi)*cos(theta)                              ];
    
    % Euler Rate Matrix
    Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

    % Body to Inertial linear velocity     
    dp = R * v; 
    % Linear velocity dynamics (Newton's 2nd law in Body Frame)
    dv = -S_omega*v + R.'*[0; 0; -g] + [0; 0; T/m]-(D*v)/m;
    % Body to Inertial  angular velocity
    dlambda = Q * omega;
    % Angular velocity dynamics (Euler's Equation)
    domega = J \ (-S_omega * J * omega + np); 
    % o Matlab recomenda fazer "J \ (-skew(w) * J * omega + np)" em vez de
    % "inv(J) * (-skew(w) * J * omega + np)" 
    
    % Group derivatives
    dX = [dp; dv; dlambda; domega];
end

function [u, v, w, theta, phi, T] = State_input_equilibrium_functions(vx, vy, psi, m, g, dx, dy, dz)
    % Shorthand for internal derivations
    kx = (m * g) ./ dx;
    ky = (m * g) ./ dy;
    
    % Longitudinal and Lateral velocity components
    V = vx .* cos(psi) + vy .* sin(psi);
    W = -vx .* sin(psi) + vy .* cos(psi);
    
    % Theta (in radians)
    % Using atan2 for numerical stability and range handling
    theta = atan2(V, kx);
    
    % u
    % Derived as u = kx * V / sqrt(kx^2 + V^2)
    u = (kx .* V) ./ sqrt(kx.^2 + V.^2);
    
    % Auxiliary term D used for phi, v, w, and T
    D = (kx .* ky + V.^2) ./ sqrt(kx.^2 + V.^2);
    
    % Phi (in radians)
    % Derived from tan(phi) = -W / D
    phi = atan2(-W, D);
    
    % Trigonometric components for final outputs
    cos_theta = kx ./ sqrt(kx.^2 + V.^2);
    cos_phi   = D ./ sqrt(W.^2 + D.^2);
    sin_phi   = -W ./ sqrt(W.^2 + D.^2);
    
    % v
    % Derived from v = -ky * sin(phi) * cos(theta)
    v = -ky .* sin_phi .* cos_theta;
    
    % w
    % Derived from w = sqrt(W^2 + D^2) - ky * cos(theta) * cos(phi)
    w = sqrt(W.^2 + D.^2) - ky .* cos_theta .* cos_phi;
    
    % T
    % Derived from T = (mg - dz * ky) * cos(theta) * cos(phi) + dz * sqrt(W^2 + D^2)
    T = (m*g - dz*ky) .* cos_theta .* cos_phi + dz .* sqrt(W.^2 + D.^2);
end

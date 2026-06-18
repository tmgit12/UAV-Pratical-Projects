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
param_values = [m_num, g_num, dx_num, dy_num, dz_num, T_op, phi_op, ...
    theta_op, psi_op];

%% 1.1 Simplified Nonlinear model
% Based on the results of Project 1, obtain a simplified nonlinear model of
% the Crazyflie drone that receives as input the vector u_lambda \in R^44,
% defined as u_lamda=[T lamda^T]^T , with the total thrust, T, and the 
% Euler angle vector lamda \in R^3, composed of the roll angle, phi, pitch 
% angle theta, and yaw angle, psi. For simplicity of the following work, 
% use consider using the linear velocity described in the inertial frame.

% (Function drone_nonlinear_dynamics is defined at the bottom)

%% 1.2 Small roll and pitch angles approximation, zero yaw
%  For small roll and pitch angles, and assuming zero yaw, obtain an 
% equivalent actuation in acceleration, u_a \in R^3, combining the effects 
% of the total thrust, gravity, and drone attitude.

D_sym = diag([dx, dy, dz]); 

% Symbolic u_lambda
u_lambda_sym = [T, phi, theta, psi];

% Symbolic rotation and equivalent acceleration
R_small_sym = small_angles_0_yaw_R(u_lambda_sym);
u_a_eq_sym  = u_lamdaTOu_a(u_lambda_sym, m, g);

%% 1.3 Linear Model
% Define the linear model for the linear motion of the drone, considering 
% the control variable u_a and the state vector x =[p^T v^T]^T.

% Symbolic dynamics
dX_sym = drone_small_angle_dynamics(X_sym, u_a_sym, m, D_sym, R_small_sym);

% Symbolic Jacobians
A_sym = jacobian(dX_sym, X_sym);
B_sym = jacobian(dX_sym, u_a_sym);

% Substitute the global numeric values into the symbolic Jacobians
A = double(subs(A_sym, param_syms, param_values));
B = double(subs(B_sym, param_syms, param_values));

%% 1.4 LQR controller
% Design an LQR controller for this linear model and test it in simulation.
% Note that you should also test the controller considering the nonlinear
% model.

% Design LQR controler
Qp = diag([120, 120, 250]);
Qv = diag([2, 2, 6]);
Q = [Qp,  zeros(3,3);
     zeros(3,3), Qv];
r = 5;
R_lqr = r * eye(3);

[K,P,poles] = lqr(A,B,Q,R_lqr);
lbds = eig(A-B*K);

%----------------------------Test models-----------------------------------
% Simulation parameters
Dt = 0.01;            % Time step
t = 0:Dt:10;          % Simulate for 10 seconds
N = length(t);

% State arrays (px, py, pz, vx, vy, vz)
x_lin = zeros(6, N+1); % Linear model state history
x_nonlin = zeros(6, N+1); % Nonlinear model state history

% Input logging arrays
u_a_log = zeros(3, N); % Log acceleration commands
u_lambda_log = zeros(4, N); % Log physical commands [T; phi; theta; psi]
x_ref_log = zeros(6, N);    % Log the moving reference target

% Numeric Drag Matrix
D_num = diag([dx_num, dy_num, dz_num]);

% Simulation "time" loop
for k = 1:N
    % 1. Define dynamic trajectory reference (Starts at t=1s)
    if k*Dt >= 1.0
        % Time elapsed since the movement started
        t_move = (k*Dt) - 1.0; 
        
        % Target Reference
        % x_ref=[5; 5; 2; 0; 0; 0]; (e.g., fly to x=5, y=5, z=2)
        x_bar_k=[5+1.5*t_move; 5; 2; 1.5; 0; 0]; % Fly at a constant 1.5 m/s
        % in the X direction. Maintain 2m altitude.
    else
        % Hovering at origin for the first 1 second
        x_bar_k = zeros(6,1);
    end

    % Log the reference so we can plot it later
    x_ref_log(:, k) = x_bar_k;
    
    % =====================================================================
    % LINEAR SYSTEM SIMULATION
    % =====================================================================
    e_x_lin = x_lin(:, k) - x_bar_k;    % State error
    u_a_lin = -K * e_x_lin;           % LQR control law
    
    % Linear dynamics: dot_x = A*x + B*u
    dot_x_lin = A * x_lin(:, k) + B * u_a_lin;
    
    % Euler integration
    x_lin(:, k+1) = x_lin(:, k) + Dt * dot_x_lin;
    
    
    % =====================================================================
    % NONLINEAR SYSTEM SIMULATION
    % =====================================================================
    e_x_nonlin = x_nonlin(:, k) - x_bar_k; % State error
    u_a_cmd    = -K * e_x_nonlin;        % LQR equivalent acceleration cmd
    
    % Transform u_a into u_lamda
    u_lambda = u_aTOu_lambda(u_a_cmd, m_num, g_num);
    
    % Nonlinear dynamics
    dot_x_nonlin = drone_nonlinear_dynamics(x_nonlin(:, k), u_lambda,...
        m_num, g_num, D_num);
    
    % Euler integration
    x_nonlin(:, k+1) = x_nonlin(:, k) + Dt * dot_x_nonlin;
    
    % Log data for plotting
    u_a_log(:, k)      = u_a_cmd;
    u_lambda_log(:, k) = u_lambda;
end

% Remove the extra last column from state arrays
x_lin(:, end)    = [];
x_nonlin(:, end) = [];

% ---------------- Plots for the state LQR --------------------------------

% Time markers
t_start = 1.0;
t_rise = t_start + 1.1;
t_settle = t_start + 3.5;

% Plotting Results
figure('Name', 'Closed-loop Trajectory', 'NumberTitle', 'off');

% Plot X Position
subplot(3,1,1);
plot(t, x_lin(1,:), 'b-', t, x_nonlin(1,:), 'r--', t, x_ref_log(1,:), ...
    'k:', 'LineWidth', 1.5);
hold on;

% Dynamic 5% bounds for X (based on the initial 5m step -> 0.25m margin)
plot(t, x_ref_log(1,:) + 0.25, 'g:', 'LineWidth', 1);
plot(t, x_ref_log(1,:) - 0.25, 'g:', 'LineWidth', 1);

xline(t_rise, 'm-.', 'Rise Time', 'LabelVerticalAlignment', 'bottom');
xline(t_settle, 'c-.', 'Settling Time', 'LabelVerticalAlignment', 'bottom');
ylabel('$$p_x$$ (m)', 'Interpreter', 'latex');
legend('Linear', 'Nonlinear', 'Target', '5% Bounds', 'Location', 'best');
title('Closed-loop Response (Position)'); grid on;

% Plot Y Position
subplot(3,1,2);
plot(t, x_lin(2,:), 'b-', t, x_nonlin(2,:), 'r--', t, x_ref_log(2,:), ...
    'k:', 'LineWidth', 1.5);
% Static bounds for Y (Target = 5m)
yline(5.25, 'g:', '+5% Overshoot');
yline(4.75, 'g:', '-5% Settling');
xline(t_rise, 'm-.', 'Rise Time', 'LabelVerticalAlignment', 'bottom');
xline(t_settle, 'c-.', 'Settling Time', 'LabelVerticalAlignment', 'bottom');
ylabel('$$p_y$$ (m)', 'Interpreter', 'latex'); grid on;

% Plot Z Position
subplot(3,1,3);
plot(t, x_lin(3,:), 'b-', t, x_nonlin(3,:), 'r--', t, x_ref_log(3,:), ...
    'k:', 'LineWidth', 1.5);
% Static bounds for Z (Target = 2m)
yline(2.1, 'g:', '+5% Overshoot');
yline(1.9, 'g:', '-5% Settling');
xline(t_rise, 'm-.', 'Rise Time', 'LabelVerticalAlignment', 'bottom');
xline(t_settle, 'c-.', 'Settling Time', 'LabelVerticalAlignment', 'bottom');
ylabel('$$p_z$$ (m)', 'Interpreter', 'latex');
xlabel('Time [s]'); grid on;

% % 3D Trajectory Plot
% figure('Name', '3D Flight Trajectory', 'NumberTitle', 'off');
% 
% % Plot the trajectories
% plot3(x_lin(1,:), x_lin(2,:), x_lin(3,:), 'b-', 'LineWidth', 1.5);
% hold on;
% plot3(x_nonlin(1,:), x_nonlin(2,:), x_nonlin(3,:), 'r--', 'LineWidth',...
%     1.5);
% plot3(x_ref_log(1,:), x_ref_log(2,:), x_ref_log(3,:), 'k:', 'LineWidth',...
%     1.5);
% % Add Start and End markers (Helps visualize the direction of flight)
% plot3(x_lin(1,1), x_lin(2,1), x_lin(3,1), 'go', 'MarkerSize', 8, ...
%     'MarkerFaceColor', 'g'); % Start Point
% plot3(x_ref_log(1,end), x_ref_log(2,end), x_ref_log(3,end), 'ro', ...
%     'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Target End Point
% % Formatting and Labels
% grid on;
% xlabel('$$p_x$$ (m)', 'Interpreter', 'latex');
% ylabel('$$p_y$$ (m)', 'Interpreter', 'latex');
% zlabel('$$p_z$$ (m)', 'Interpreter', 'latex');
% title('3D Trajectory Tracking (Linear vs Nonlinear)');
% legend('Linear', 'Nonlinear', 'Target Path', 'Start', 'Final Target', ...
%     'Location', 'best');
% % Set the viewing angle
% view(-45, 30); % Sets a nice isometric viewing angle (Azimuth, Elevation)
% axis equal;


% ---------------- Plots for LQR Actuation --------------------------------
figure('Name', 'Nonlinear Actuation', 'NumberTitle', 'off');

% Plot Thrust
subplot(2,1,1);
plot(t, u_lambda_log(1,:), 'g', 'LineWidth', 1.5);
ylabel('Thrust (N)');
title('Closed-loop Actuation Commands'); grid on;

% Plot Roll and Pitch Angles
subplot(2,1,2);
plot(t, rad2deg(u_lambda_log(2,:)), 'b', t, ...
   rad2deg(u_lambda_log(3,:)), 'r', 'LineWidth', 1.5);
ylabel('Angle (deg)');
legend('\phi (Roll)', '\theta (Pitch)', 'Location', 'best');
xlabel('Time [s]'); grid on;

%% 1.5 Consider now an error state vector defined as ˜x = x - _x, where we 
% assume that the reference state, _x, is driven by the same dynamics as x.
% Obtain the equivalent state space model.

A_err = A;
B_err = B;
%% 1.6 Design an LQR controller for this error model and test it in 
% simulation. Comment on the difference between this linear model and the 
% previous one.

% Design error LQR controler
Qp_err = diag([350, 350, 400]);
Qv_err = diag([30, 30, 35]);    
Q_err = [Qp_err, zeros(3,3);
    zeros(3,3), Qv_err];
 r_err = 5;
R_lqr_err = r_err * eye(3);

% Recalculate LQR gain
[K_err, P_err, poles_err] = lqr(A_err, B_err, Q_err, R_lqr_err);

% Closed-loop error dynamics: x_tilde_dot = (A - B*K) * x_tilde
A_cl_err = A_err - B_err * K_err;

% State-space
C_err = [eye(3), zeros(3,3)]; % observe error state
D_err = zeros(3,3);           % No direct feedthrough
sys_cl_error = ss(A_cl_err, B_err, C_err, D_err);

%----------------------------Test models-----------------------------------


x_tilde = zeros(6, N); % Error state history
x_bar = zeros(6, N); % Reference state history
x = zeros(6, N+1); % Actual state history

u_a_tilde_log = zeros(3, N); % Control input history

x(:,1) = zeros(6,1); % Initial condition

x_tilde(:,1) = x(:,1)-x_bar(:,1); % x_tilde = x - x_bar

for k = 1:N
    % Reference trajectory x_bar
    if k*Dt >= 1.0
        t_move = (k*Dt) - 1.0;
        % Target Reference (e.g., fly to x=5, y=5, z=2)
        % x_ref=[5; 5; 2; 0; 0; 0]; 

        % Fly at a constant 1.5 m/s in the X direction. Maintain 2m 
        % altitude.
        x_bar_k=[5+1.5*t_move; 5; 2; 1.5; 0; 0]; 
    else
        % Hover at origin
        x_bar_k = zeros(6,1);
    end

    % Log reference
    x_bar(:,k) = x_bar_k;

    % Calculat error state
    x_tilde(:,k) = x(:,k) - x_bar_k;

    % LQR control law on the error model
    u_a_tilde = -K_err * x_tilde(:,k);

    % Error dynamics
    % x_tilde_dot = (A - B*K)x_tilde
    x_tilde_dot = A * x_tilde(:, k) + B * u_a_tilde;
    %x_tilde_dot = A_cl_err * x_tilde(:,k);

    % Integrate error dynamics
    x_tilde(:,k+1) = x_tilde(:,k) + Dt * x_tilde_dot;

    % Recover actual state from x = x_tilde + x_bar
    x(:,k+1) = x_tilde(:,k+1) + x_bar_k;

    % Log control input
    u_a_tilde_log(:,k) = u_a_tilde;
end

% Remove extra sample
 x(:,end) = [];
 x_tilde(:,end) = [];

% ---------------- Plots for the Error-Model LQR --------------------------

figure('Name','Closed-loop Tracking Comparison','NumberTitle','off');

% Plot X Position
subplot(3,1,1);
plot(t, x_lin(1,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(1,:), 'r--', 'LineWidth',1.5);
plot(t, x(1,:), 'g-.', 'LineWidth',1.5);
plot(t, x_bar(1,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_x$$ (m)', 'Interpreter','latex');
legend('Linear Model', ...
       'Nonlinear Model', ...
       'Error-State Linear Model', ...
       'Reference Trajectory', ...
       'Location','best');
title('Trajectory Tracking Comparison');
grid on;

% Plot Y Position
subplot(3,1,2);
plot(t, x_lin(2,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(2,:), 'r--', 'LineWidth',1.5);
plot(t, x(2,:), 'g-.', 'LineWidth',1.5);
plot(t, x_bar(2,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_y$$ (m)', 'Interpreter','latex');
grid on;

% Plot Z Position
subplot(3,1,3);
plot(t, x_lin(3,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(3,:), 'r--', 'LineWidth',1.5);
plot(t, x(3,:), 'g-.', 'LineWidth',1.5);
plot(t, x_bar(3,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_z$$ (m)', 'Interpreter','latex');
xlabel('Time [s]');
grid on;


% Error Plot

figure('Name','Error-State Evolution','NumberTitle','off');

subplot(3,1,1);
plot(t, x_tilde(1,:), 'LineWidth',1.5);
ylabel('$$\tilde{p}_x$$ (m)','Interpreter','latex');
title('Error-State Convergence');
grid on;

subplot(3,1,2);
plot(t, x_tilde(2,:), 'LineWidth',1.5);
ylabel('$$\tilde{p}_y$$ (m)','Interpreter','latex');
grid on;

subplot(3,1,3);
plot(t, x_tilde(3,:), 'LineWidth',1.5);
ylabel('$$\tilde{p}_z$$ (m)','Interpreter','latex');
xlabel('Time [s]');
grid on;


% 3D Tracking Plot

figure('Name','3D Error-Model Tracking','NumberTitle','off');

plot3(x(1,:), x(2,:), x(3,:), 'b', 'LineWidth',1.5);
hold on;
plot3(x_bar(1,:), x_bar(2,:), x_bar(3,:), 'k:', 'LineWidth',1.5);
% Start point
plot3(x(1,1), x(2,1), x(3,1), 'go', 'MarkerSize',8, 'MarkerFaceColor','g');
% Final reference point
plot3(x_bar(1,end), x_bar(2,end), x_bar(3,end), 'ro', 'MarkerSize',8, ...
    'MarkerFaceColor','r');
grid on;
axis equal;
xlabel('$$p_x$$ (m)','Interpreter','latex');
ylabel('$$p_y$$ (m)','Interpreter','latex');
zlabel('$$p_z$$ (m)','Interpreter','latex');
title('3D Trajectory Tracking with Error-State LQR');
legend('Actual Trajectory $$x$$', 'Reference Trajectory $$\bar{x}$$', ...
       'Start', 'Final Reference', 'Interpreter','latex', 'Location', ...
       'best');
view(-45,30);

%% 2.1 Design a nonlinear controller with actuation in body accelerations, 
% u_a \in R^3 and assuming a zero yaw, that is able to achieve asymptotic 
% stability in the Lyapunov sense (prove this result).

rise_time = 1.1;
setlling_time = 3.2;
overshoot = 5;

fator_amortecimento = sqrt(log(overshoot/100)^2/(pi^2+log(overshoot/100)^2));
freq_natural_risetime = (3.7*fator_amortecimento) / rise_time; %antes era 3.7 mas queria deixar 4x mais rápido
freq_natural_setlling_time = 3/(setlling_time*fator_amortecimento); %antes era 3 mas queria deixar 1.48x mais rápido

freq_natural = max(freq_natural_setlling_time,freq_natural_risetime);

Kp_gain = 2 * fator_amortecimento * freq_natural; %o fator de amortecimento é 5.30 no eixo Z
Kd_gain = freq_natural^2; % tem de ser 2.83 no eixo Z

%Kp_mat = diag([Kp_gain, Kp_gain, Kp_gain]); % Position error gains
%Kv_mat = diag([Kd_gain, Kd_gain, Kd_gain]);

Kp_mat = diag([16, 16.0, 30.0]); % Z ligeiramente mais forte
Kv_mat = diag([12, 12, 8]);

%% 2.2 Test this controller in simulation and compare with the previous 
% linear controller.

x_lyap = zeros(6, N); % Nonlinear model state history
x_tilde = zeros(6, N); % Error state history x_tilde = x - x_bar
x_bar = zeros(6, N); % Reference state history


% Input logging arrays
u_a_log = zeros(3, N); % Log acceleration commands
u_lambda_log = zeros(4, N); % Log physical commands [T; phi; theta; psi]
x_ref_log = zeros(6, N);    % Log the moving reference target

for k = 1:N
    % Reference trajectory x_bar
    if k*Dt >= 1.0
        t_move = (k*Dt) - 1.0;
        % Target Reference
        % x_ref=[5; 5; 2; 0; 0; 0]; (e.g., fly to x=5, y=5, z=2)
        % Fly at a constant 1.5 m/s in the X direction. Maintain 2m 
        % altitude.
        x_bar_k=[5+1.5*t_move; 5; 2; 1.5; 0; 0]; 
    else
        % Hover at origin
        x_bar_k = zeros(6,1);
    end

    % Log reference
    x_bar(:,k) = x_bar_k;

    % Calculat error state
    x_tilde(:,k) = x_lyap(:,k) - x_bar_k;
    e_p = x_tilde(1:3, k); % Position error
    e_v = x_tilde(4:6, k); % Velocity error

    v_current = x_lyap(4:6, k);

    if k == 1
        R_current = eye(3);
    else
        % Z-Y-X Euler Rotation from previous u_lambda
        phi=u_lambda(2);
        theta=u_lambda(3);
        psi=0;
        R_current = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
         -sin(theta) , sin(phi)*cos(theta) , cos(phi)*cos(theta)];        
    end

    dot_v_current = [0; 0; 0];

    % NONLINEAR CONTROL LAW
    % u_a = a_ref + (R*D*R'*v)/m - Kp*e_p - Kv*e_v 
    u = - (Kp_mat * e_p) - (Kv_mat * e_v);
    u_a_cmd = dot_v_current + (R_current * D_num * R_current' * v_current)/m_num + u ;

    u_lambda = u_aTOu_lambda(u_a_cmd, m_num, g_num);
    dot_x_nonlin = drone_nonlinear_dynamics(x_lyap(:, k), u_lambda, m_num, g_num, D_num);
    x_lyap(:, k+1) = x_lyap(:, k) + Dt * dot_x_nonlin;
end

% Remove extra sample
%x(:,end) = [];
x_lyap(:,end) = [];

% ---------------- Plots for the Model Comparisons ------------------------
figure('Name','Closed-loop Tracking Comparison','NumberTitle','off');

% Plot X Position
subplot(3,1,1);
plot(t, x_lin(1,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(1,:), 'r--', 'LineWidth',1.5);
plot(t, x(1,:), 'g-.', 'LineWidth',1.5);
plot(t, x_lyap(1,:), 'm-', 'LineWidth',1.5); % <--- NEW: Lyapunov Model
plot(t, x_bar(1,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_x$$ (m)', 'Interpreter','latex');
% Updated Legend
legend('Linear Model (LQR)', ...
       'Nonlinear Model (LQR)', ...
       'Error-State Linear (LQR)', ...
       'Nonlinear Model (Lyapunov)', ...
       'Reference Trajectory', ...
       'Location','best');
title('Trajectory Tracking Comparison');
grid on;

% Plot Y Position
subplot(3,1,2);
plot(t, x_lin(2,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(2,:), 'r--', 'LineWidth',1.5);
plot(t, x(2,:), 'g-.', 'LineWidth',1.5);
plot(t, x_lyap(2,:), 'm-', 'LineWidth',1.5); % <--- NEW: Lyapunov Model
plot(t, x_bar(2,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_y$$ (m)', 'Interpreter','latex');
grid on;

% Plot Z Position
subplot(3,1,3);
plot(t, x_lin(3,:), 'b-',  'LineWidth',1.5); hold on;
plot(t, x_nonlin(3,:), 'r--', 'LineWidth',1.5);
plot(t, x(3,:), 'g-.', 'LineWidth',1.5);
plot(t, x_lyap(3,:), 'm-', 'LineWidth',1.5); % <--- NEW: Lyapunov Model
plot(t, x_bar(3,:), 'k:',  'LineWidth',1.5);
ylabel('$$p_z$$ (m)', 'Interpreter','latex');
xlabel('Time [s]');
grid on;

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
         -sin(theta) , sin(phi)*cos(theta) , cos(phi)*cos(theta)];
         
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

%--------------------------------------------------------------------------
% For 1.4 LRQ control
function u_lambda = u_aTOu_lambda(u_a, m, g)
    uax = u_a(1); uay = u_a(2); uaz = u_a(3);
    
    % Total thrust (compensating gravity + desired Z acceleration)
    T = m * (uaz + g);
    if T < 0.01; T = 0.01; end % Prevent negative thrust
    
    % Pitch and Roll commands based on LQR X/Y acceleration
    theta =  (m * uax) / T;
    phi   = -(m * uay) / T;
    
    % Angle Saturation (drone will crash if requested to tilt 90 deg)
    max_angle = 30 * (pi/180); % 30 degrees limit
    phi   = max(min(phi, max_angle), -max_angle);
    theta = max(min(theta, max_angle), -max_angle);
    
    u_lambda = [T; phi; theta; 0]; % Yaw is zero
end

%--------------------------------------------------------------------------
% For 2.1
function u_a_nonlinear = u_a_nonlinear(u_a, m, g)
    
end

% =========================================================================
% Crazyflie 2.1 Nonlinear Model Simulation
% =========================================================================
clear; clc; close all;

% Crazyflie 2.1 Parameters 
m = 0.027; % Mass in kg
g = 9.81;  % Gravity in m/s^2
l = 0.05;  % Arm length in m
Jx = 3.3e-5; Jy = 3.6e-5; Jz = 5.9e-5; 
J = diag([Jx, Jy, Jz]); % Inertia Matrix in kg*m^2
Fax = 0.01; Fay = 0.01; Faz = 0.01;
D = diag([Fax, Fay, Faz]); % Drag Matrix in kg/s

%% Simulation
dt=0.01;
t = 0:dt:5; % Simulation interval
X = zeros(12, length(t));
u = zeros(4,length(t));  % u = [T; n_p] (4 inputs: thrust + 3 moments)

% Initial state vector: [p; v; lambda; omega]
% p=[14;15;16], v=[0;0;0], lambda=[0;0;0], omega=[0;0;0]  
X(1:3, 1) = [14; 15; 16];

% Inputs: Small constant thrust, equal to all rotors
% T_Total = T1+T2+T3+T4
T_total = m * g; 
n_p = [0; 0; 0]; % Equal thrust means zero moments
u = [T_total; n_p] * ones(1, length(t));



for i = 1:(length(t) - 1)
    current_t = t(i);
    current_X = X(:, i);
    current_u = u(:, i);

    dX = drone_dynamics(current_X, m, g, J, current_u, D);

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
legend('v_x', 'v_y', 'v_z');

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
    
    % Rotation Matrix (Z-Y-X Euler to Body->Inertial)
    R = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
         cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
         -sin(theta)        , sin(phi)*cos(theta)                             , cos(phi)*cos(theta)                              ];
    
    % Euler Rate Matrix
    Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
         
    dp = R * v; % Body to Inertial linear velocity
    dv = -skew(omega)*v + R.'*[0; 0; -g] + [0; 0; T/m]-(D*v)/m; % Linear velocity dynamics (Newton's 2nd law in Body Frame)
    dlambda = Q * omega; % Body to Inertial  angular velocity
    domega = J \ (-skew(omega) * J * omega + np); % Angular velocity dynamics (Euler's Equation)
    % o Matlab recomenda fazer "J \ (-skew(w) * J * omega + np)" em vez de
    % "inv(J) * (-skew(w) * J * omega + np)" 
    
    % Group derivatives
    dX = [dp; dv; dlambda; domega];
end

clear; clc; close all;
    
syms px py pz real % derivada da posicao
syms vx vy vz real % derivada da velocidade
syms phi theta psi real % angulos
syms wx wy wz real % velocidade angular
states_sym = [px; py; pz; ...
              vx; vy; vz; ...  
              phi; theta; psi; ... 
              wx; wy; wz];

syms Jx Jy Jz real % parametros da matriz de Inercia
syms m g fax fay faz real % parametros de massa, gravidade e atrito

syms T npx npy npz  real % input

J = [Jx, 0, 0;
     0, Jy, 0;
     0, 0, Jz];

% Rotation Matrix (Z-Y-X Euler to Body->Inertial)
Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
R = Rz * Ry * Rx; % Body to Inertial

% Euler Rate Matrix
Q = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
    0, cos(phi),           -sin(phi);
    0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

% Vetores
p = [px; py; pz];
v = [vx; vy; vz];
lambda = [phi; theta; psi];
omega = [wx; wy; wz];

x=[p; v; lambda; omega];

np = [npx; npy; npz];
fg_inertial = [0; 0; -m*g];

D = [fax, 0, 0;
     0, fay, 0;
     0, 0, faz];


% f(x,u)
dp = R * v;
dv = cross(-omega, v) + (1/m) * (R' * fg_inertial - D*v + [0; 0; T]);
dlambda = Q * omega;
domega = J \ (cross(-omega, J * omega) + sym(np));



%% Matrix A
% A=df/dx ; x=[p; v; lambda; omega]

dp_p=jacobian(dp,p);
dp_v=jacobian(dp,v);
dp_lambda=jacobian(dp,lambda);
dp_omega=jacobian(dp,omega);

dv_p=jacobian(dv,p);
dv_v=jacobian(dv,v);
dv_lambda=jacobian(dv,lambda);
dv_omega=jacobian(dv,omega);

dlambda_p = jacobian(dlambda, p);
dlambda_v = jacobian(dlambda, v);
dlambda_lambda = jacobian(dlambda,lambda);
dlambda_omega = jacobian(dlambda,omega);

domega_p = jacobian(domega, p);
domega_v = jacobian(domega, v);
domega_lambda = jacobian(domega, lambda);
domega_omega = jacobian(domega, omega);

A = [dp_p, dp_v, dp_lambda, dp_omega;
     dv_p, dv_v, dv_lambda, dv_omega;
     dlambda_p, dlambda_v, dlambda_lambda, dlambda_omega;
     domega_p, domega_v, domega_lambda, domega_omega];

%% Matrix B
%B=df/du 
u = [T; np];

dp_u=jacobian(dp,u);
dv_u=jacobian(dv,u);
dlambda_u = jacobian(dlambda, u);
domega_u = jacobian(domega, u);

B = [dp_u; dv_u; dlambda_u; domega_u];

%% Linearização
params_sym = [m; g; Jx; Jy; Jz; fax; fay; faz];
params_val = [0.027; 9.81; 3.3e-5; 3.6e-5; 5.9e-5; 0.01; 0.01; 0.01];

all_syms = [states_sym; params_sym; T; npx; npy; npz];

%% OP1
% O mais logico em hover e phi_e=0 e theta_e=0, psi_e pode ser qualquer
% valor, mexer nos valores.
xe_1=[10;10;10; 0;0;0; 0;0;0; 0;0;0]; 
ue_1=[params_val(1)*params_val(2)/(cos(xe_1(7))*cos(xe_1(8)));0;0;0];

A_OP1 = double(subs(A, all_syms, [xe_1; params_val; ue_1]));
B_OP1 = double(subs(B, all_syms, [xe_1; params_val; ue_1]));

x_til1=x-xe_1;
u_til1=u-ue_1;

dx_til_OP1=A_OP1*x_til1+B_OP1*u_til1;

%% A_OP2
% Primeiro simplificar com drone a andar so no eixo do x, vy=0 e phi=0
theta_e2 = pi/4;
phi_e2   = 0;

vx_e2 = (params_val(1) * params_val(2) * sin(pi/4)) / params_val(6);
vy_e2 = (-params_val(1) * params_val(2) * sin(phi_e2) * cos(theta_e2))...
    / params_val(7);
vz_e2 = 0;

xe_2 = [10; 10; 10; vx_e2; vy_e2; vz_e2; phi_e2; theta_e2; 0; 0; 0; 0];
ue_2 = [params_val(1) * params_val(2) * cos(phi_e2) * cos(theta_e2); ...
    0; 0; 0];

A_OP2 = double(subs(A, all_syms, [xe_2; params_val; ue_2]));
B_OP2 = double(subs(B, all_syms, [xe_2; params_val; ue_2]));

x_til2=x-xe_2;
u_til2=u-ue_2;
dx_til_OP2=A_OP2*x_til2+B_OP2*u_til2;

%% Plot Eigenvalues

figure;
hold on;
grid on;

e1 = eig(A_OP1);
e2 = eig(A_OP2);

h1 = plot(real(e1), imag(e1), 'ro', 'MarkerSize', 10, 'LineWidth', ...
    1.5, 'DisplayName', 'OP1 (Hover)');
h2 = plot(real(e2), imag(e2), 'bx', 'MarkerSize', 10, 'LineWidth', ...
    1.5, 'DisplayName', 'OP2 (Horizontal)');

xline(0, 'k-', 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off');
xlabel('Re');
ylabel('Imag');
legend([h1, h2], 'Location', 'best');

%% Models controlability, observability and stability
n_OP1 = length(A_OP1);
n_OP2 = length(A_OP2);
C_OP1 = eye(12);
C_OP2 = eye(12);

% OP1
%   controlability
Mc_OP1 = ctrb( A_OP1, B_OP1);  % Mc = [ B | A*B | ... | A^(n-1)*B ]
rank(Mc_OP1)
if ( rank(Mc_OP1) == n_OP1 )
    disp( 'OP1 linear model is fully controllable' )
else
    disp( 'OP1 linear model is not fully controllable')
end

%   observability
Mo_OP1 = obsv(A_OP1 , C_OP1);  % Mo = [ C; C*A ]

rank(Mo_OP1)
if ( rank(Mo_OP1) == n_OP1 )
    disp( 'OP1 linear model is fully observable' )
else
    disp( 'OP1 linear model is not fully observable')
end

%   stability
J_OP1 = jordan(A_OP1);

% OP2
%   controlability
Mc_OP2 = ctrb( A_OP2, B_OP2);  % Mc = [ B | A*B | ... | A^(n-1)*B ]
rank(Mc_OP2)
if ( rank(Mc_OP2) == n_OP2 )
    disp( 'OP2 linear model is fully controllable' )
else
    disp( 'OP2 linear model is not fully controllable')
end

%   observability
Mo_OP2 = obsv(A_OP2 , C_OP2);  % Mo = [ C; C*A ]

rank(Mo_OP2)
if ( rank(Mo_OP2) == n_OP2 )
    disp( 'OP2 linear model is fully observable' )
else
    disp( 'OP2 linear model is not fully observable')
end

%   stability
J_OP2 = jordan(A_OP2);

%% Transfer functions P13
% OP1

sys_OP1 = ss(A_OP1, B_OP1, eye(12), zeros(12,4));
sys_tf_OP1 = tf(sys_OP1);

G_nx_phi_sym_OP1 = sys_tf_OP1(7, 2); % nx/phi
G_ny_theta_sym_OP1 = sys_tf_OP1(8, 3); % ny/theta
G_nz_psi_sym_OP1 = sys_tf_OP1(9, 4); % nz/psi

G_px_theta_sym_OP1 = sys_tf_OP1(1, 3) / sys_tf_OP1(8, 3); % px/theta
G_py_phi_sym_OP1 = sys_tf_OP1(2, 2) / sys_tf_OP1(7, 2); % py/phi
G_pz_T_sym_OP1 = sys_tf_OP1(3, 1); % pz/T

figure;
rlocus(G_pz_T_sym_OP1);
figure;
rlocus(G_py_phi_sym_OP1);


% OP2
sys_OP2 = ss(A_OP2, B_OP2, eye(12), zeros(12,4));
sys_tf_OP2 = tf(sys_OP2);

G_nx_phi_sym_OP2 = sys_tf_OP2(7, 2); % nx/phi
G_ny_theta_sym_OP2 = sys_tf_OP2(8, 3); % ny/theta
G_nz_psi_sym_OP2 = sys_tf_OP2(9, 4); % nz/psi

G_px_theta_sym_OP2 = sys_tf_OP2(1, 3) / sys_tf_OP2(8, 3); % px/theta
G_py_phi_sym_OP2 = sys_tf_OP2(2, 2) / sys_tf_OP2(7, 2); % py/phi
G_pz_T_sym_OP2 = sys_tf_OP2(3, 1); % pz/T

figure;
rlocus(G_pz_T_sym_OP2);
figure;
rlocus(G_py_phi_sym_OP2);

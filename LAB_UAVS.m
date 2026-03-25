clear; clc; close all;
    
syms px py pz real % derivada da posição
syms vx vy vz real % derivada da velocidade
syms phi theta psi real % ângulos
syms wx wy wz real % velocidade angular
states_sym = [px; py; pz; ...
              vx; vy; vz; ...  
              phi; theta; psi; ... 
              wx; wy; wz];

syms Jx Jy Jz real % parâmetros da matriz de Inércia
syms m g fax fay faz real % parâmetros de massa, gravidade e atrito

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

% matrizes
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
dv = cross(-omega, v) + (1/m) * (R' * fg_inertial + D*v + [0; 0; T]);
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


%% A_OP1
% O mais lógico em hover é phi_e=0 e theta_e=0, psi_e pode ser qualquer
% valor, mexer nos valores.
params_sym_1 = [m; g; Jx; Jy; Jz; fax; fay; faz];
params_val_1 = [0.027; 9.81; 3.3e-5; 3.6e-5; 5.9e-5; 0.01; 0.01; 0.01];

input_val_1 = [0; 0; 0; 0];
input_sym_1 = [T; npx; npy; npz];

all_syms = [states_sym; params_sym_1; input_sym_1];

xe_1=[10;10;10; 0;0;0; 0;pi/4;0; 0;0;0]; 
ue_1=[-m*g*cos(xe_1(7))*cos(xe_1(8));0;0;0];
all_vals_1 = [xe_1; params_val_1; input_val_1];

A_OP1 = double(subs(A, all_syms, all_vals_1));
B_OP1 = double(subs(B, all_syms, all_vals_1));
eig_OP1=eig(A_OP1);

figure;
plot(real(eig_OP1), imag(eig_OP1), 'ro', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('OP1');

x_til1=x-xe_1;
u_til1=u-ue_1;

dx_til_OP1=A_OP1*x_til1+B_OP1*u_til1;

%% A_OP2
params_sym_2 = [m; g; Jx; Jy; Jz; fax; fay; faz];
params_val_2 = [0.027; 9.81; 3.3e-5; 3.6e-5; 5.9e-5; 0.01; 0.01; 0.01];

input_sym_2 = [T; npx; npy; npz];
input_val_2 = [0; 0; 0; 0];

all_syms_2 = [states_sym; params_sym_2; input_sym_2];

xe_2=[10;10;10; vx;vy;0; 0;pi/4;0; 0;0;0];

vx_2 = m * g * sin(xe_2(8)) / fax;
vy_2 = (-m * g * cos(xe_2(8)) * cos(xe_2(7))) / fay;
vx_2 = subs(vx_2,params_sym_2, params_val_2);
vy_2 = subs(vy_2,params_sym_2, params_val_2);
xe_2=subs(xe_2, [vx;vy], [vx_2; vy_2]);

ue_2=[m*g*cos(xe_2(7))*cos(xe_2(8));0;0;0];
ue_2=subs(ue_2, params_sym_2, params_val_2);
all_vals_2 = [xe_2; params_val_2; ue_2];


A_OP2 = double(subs(A, all_syms_2, all_vals_2));
B_OP2 = double(subs(B, all_syms_2, all_vals_2));
eig_OP2=eig(A_OP2);

figure;
plot(real(eig_OP2), imag(eig_OP2), 'ro', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('OP2');


x_til2=x-xe_2;
u_til2=u-ue_2;
dx_til_OP2=A_OP2*x_til2+B_OP2*u_til2;
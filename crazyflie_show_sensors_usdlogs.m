% crazyflie load usd card log csv file for sensor analysis
clear all

% available files:
csvfilename_motors_off = 'L2Data1.csv';
csvfilename_hover = 'L2Data2.csv'; % OP1
csvfilename_motion_x = 'L2Data3.csv'; % OP2
csvfilename_motion_y = 'L2Data4.csv';
csvfilename_motion_z = 'L2Data5.csv';
csvfilename_motion_inf = 'L2Data6.csv';

%% OP1
% read file
csvfilename = csvfilename_hover;
array = dlmread(csvfilename,',',1,0);
%T = table2array(readtable(csvfilename)); % Matlab only

% get data from table
time = array(:,1)'*1e-3;
pos = array(:,2:4)'; % [m] groundtruth
vel = array(:,5:7)'; % [m/s] groundtruth
lbd = array(:,8:10)'; % [deg] groundtruth
gyro = array(:,11:13)'; % [deg/s] sensor
acc = array(:,14:16)'; % [Gs] sensor
baro_asl = array(:,17)'; % [m] sensor
% lighthouse = array(:,18:20))'; % [m]

% convert date to print format
t = time - time(1);

% plot data
initPlots;
crazyflie_show_sensor_data(t,acc,gyro,baro_asl,pos,vel,lbd);

% 2.1

% X axis

px = pos(1,:);
vx = vel(1,:);
phi = lbd(1,:);
wx = gyro(1,:);
ax = acc(1,:);

% Y axis

py = pos(2,:);
vy = vel(2,:);
theta = lbd(2,:);
wy = gyro(2,:);
ay = acc(2,:);

% Z axis

pz = pos(3,:);
vz = vel(3,:);
psi = lbd(3,:);
wz = gyro(3,:);
az = acc(3,:);

% Modos
Ts = mean(diff(t)); % Sampling time from your time vector

% Input: theta, Output: px

data_x = iddata(px', theta', Ts);
data_x = detrend(data_x);

data_y = iddata(py', phi', Ts);
data_y = detrend(data_y);

az_T = (az - 1) * 9.81; 
data_z = iddata(pz', az_T', Ts);
data_z = detrend(data_z);


% Estimate transfer functions
np = 2; nz = 0; % nº de polos e nº de zeros
sys_x = tfest(data_x, np, nz);
sys_y = tfest(data_y, np, nz);
sys_z = tfest(data_z, np, nz); 

disp('G_px_theta:');
tf(sys_x)

disp('G_py_phi:');
tf(sys_y)

%% OP2
% read file
csvfilename = csvfilename_motion_x;
array = dlmread(csvfilename,',',1,0);
%T = table2array(readtable(csvfilename)); % Matlab only

% get data from table
time = array(:,1)'*1e-3;
pos = array(:,2:4)'; % [m] groundtruth
vel = array(:,5:7)'; % [m/s] groundtruth
lbd = array(:,8:10)'; % [deg] groundtruth
gyro = array(:,11:13)'; % [deg/s] sensor
acc = array(:,14:16)'; % [Gs] sensor
baro_asl = array(:,17)'; % [m] sensor
% lighthouse = array(:,18:20))'; % [m]

% convert date to print format
t = time - time(1);

% plot data
initPlots;
crazyflie_show_sensor_data(t,acc,gyro,baro_asl,pos,vel,lbd);

% 2.1

% X axis

px = pos(1,:);
vx = vel(1,:);
phi = lbd(1,:);
wx = gyro(1,:);
ax = acc(1,:);

% Y axis

py = pos(2,:);
vy = vel(2,:);
theta = lbd(2,:);
wy = gyro(2,:);
ay = acc(2,:);

% Z axis

pz = pos(3,:);
vz = vel(3,:);
psi = lbd(3,:);
wz = gyro(3,:);
az = acc(3,:);

% Modos
Ts = mean(diff(t)); % Sampling time from your time vector

% Input: theta, Output: px

data_x = iddata(px', theta', Ts);
data_x = detrend(data_x);

data_y = iddata(py', phi', Ts);
data_y = detrend(data_y);

az_T = (az - 1) * 9.81; 
data_z = iddata(pz', az_T', Ts);
data_z = detrend(data_z);


% Estimate transfer functions
np = 2; nz = 0; % nº de polos e nº de zeros
sys_x = tfest(data_x, np, nz);
sys_y = tfest(data_y, np, nz);
sys_z = tfest(data_z, np, nz); 


disp('G_px_theta:');
tf(sys_x)

disp('G_py_phi:');
tf(sys_y)

disp('G_pz_T:');
tf(sys_z)

%% Comparison

figure;
compare(data_x, sys_x);
title('\theta to p_x');

figure;
compare(data_y, sys_y);
title('\phi to p_y');

figure;
compare(data_z, sys_z);
title('T to p_z)');


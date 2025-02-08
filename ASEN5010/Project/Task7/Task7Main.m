%% ASEN 5010 Task 7 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all 

%% Setup
% Include all of our lovely task functions
addpath('..\..\Utilities')
addpath('..\Task1')
addpath('..\Task2')
addpath('..\Task3')
addpath('..\Task4')
addpath('..\Task5')
addpath('..\Task6')

% Make Mars
R_Mars = 3396.19; % km
[marsX, marsY, marsZ] = sphere(100);
marsX = R_Mars*marsX;
marsY = R_Mars*marsY;
marsZ = R_Mars*marsZ;

% Initial orbit parameters
h_LMO = 400; % km
h_GMO = 17028.01; % km
radius_LMO = R_Mars + h_LMO; % km
radius_GMO = R_Mars + h_GMO; % km

w_0_LMO = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_0_LMO = deg2rad([20; 30; 60]); % Omega, i, theta

w_0_GMO = [0; 0; 0.0000709003]; % rad/s, O frame coords
EA_0_GMO = deg2rad([0; 0; 250]); % Omega, i, theta

x_0_LMO = [radius_LMO; EA_0_LMO; w_0_LMO];
x_0_GMO = [radius_GMO; EA_0_GMO; w_0_GMO];

% Initial attitude parameters
sigBN_0 = [0.3; -0.4; 0.5]; % unitless
omegBN_0 = deg2rad([1.00; 1.75; -2.20]); % Converted to rad/s from deg/s
u_0 = zeros(3,1); % Nm
I = diag([10, 5, 7.5]); % kgm^2, body coords

x_0_att = {I; sigBN_0; omegBN_0; u_0};

% RK4 params
t0 = 0;
dt = 1;
tf = 1000;

%% Propagate orbits and attitude with 0 control
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

out_Attitude = RK4_Attitude(x_0_att, t0, dt, tf, 'no torque');

%% Plot simulation outputs
t = out_Attitude(:,1);
sig = out_Attitude(:,2:4);
w = out_Attitude(:,5:7);
u = out_Attitude(:,8:10);

% Orbit
figOrbit = figure;

ax2 = axes();
marsStars = imread("marsStars.jpg");
imshow(marsStars, 'parent', ax2)

ax1 = axes();
title("Orbit Simulation", 'Color', 'w')
hold on
grid on
axis equal
plot3(out_LMO(:,2), out_LMO(:,3), out_LMO(:,4), 'LineWidth', 3)
plot3(out_GMO(:,2), out_GMO(:,3), out_GMO(:,4), 'LineWidth', 3)

mars = surf(marsX, marsY, -marsZ); % Need to flip the sphere for image to map properly
marsSurface = imread("marsSurface.jpg");
set(mars,'FaceColor','texturemap','cdata',marsSurface,'edgecolor','none');

set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w')

xlabel("$\hat{n}_1$", "Interpreter","latex")
ylabel("$\hat{n}_2$", "Interpreter","latex")
zlabel("$\hat{n}_3$", "Interpreter","latex")

view([30 35])

% Attitude
figAttitude = figure;

sgtitle("Nano-satellite Attitude Evolution Over Time")

subplot(3,3,1)
hold on
grid on
title("\sigma_1 vs. time")
plot(t, sig(:,1))
xlabel("Time [sec]")
ylabel("\sigma_1")

subplot(3,3,2)
hold on
grid on
title("\omega_1 vs. time")
plot(t, w(:,1))
xlabel("Time [sec]")
ylabel("\omega_1")

subplot(3,3,3)
hold on
grid on
title("u_1 vs. time")
plot(t, u(:,1))
xlabel("Time [sec]")
ylabel("u_1")

subplot(3,3,4)
hold on
grid on
title("\sigma_2 vs. time")
plot(t, sig(:,2))
xlabel("Time [sec]")
ylabel("\sigma_2")

subplot(3,3,5)
hold on
grid on
title("\omega_2 vs. time")
plot(t, w(:,1))
xlabel("Time [sec]")
ylabel("\omega_2")

subplot(3,3,6)
hold on
grid on
title("u_2 vs. time")
plot(t, u(:,2))
xlabel("Time [sec]")
ylabel("u_2")

subplot(3,3,7)
hold on
grid on
title("\sigma_3 vs. time")
plot(t, sig(:,3))
xlabel("Time [sec]")
ylabel("\sigma_3")

subplot(3,3,8)
hold on
grid on
title("\omega_3 vs. time")
plot(t, w(:,3))
xlabel("Time [sec]")
ylabel("\omega_3")

subplot(3,3,9)
hold on
grid on
title("u_3 vs. time")
plot(t, u(:,3))
xlabel("Time [sec]")
ylabel("u_3")

%% Extract answers to text files
t = 500;

sig_ans = sig(out_Attitude(:,1) == t, :)';
w_ans = w(out_Attitude(:,1) == t, :)';

H_ans = I*w_ans; % In body coords
T_ans = 0.5*w_ans'*I*w_ans;

BN = MRP2DCM(sig_ans);
NB = BN';

H_ans_2 = NB*H_ans; % In inertial coords

f1 = fopen("H_ans.txt", "w");
ans_H = fprintf(f1, "%.3f %.3f %.3f", H_ans(1), H_ans(2), H_ans(3));
fclose(f1);

f2 = fopen("T_ans.txt", "w");
ans_T = fprintf(f2, "%.5f", T_ans);
fclose(f2);

f3 = fopen("sig_ans.txt", "w");
ans_sig = fprintf(f3, "%.3f %.3f %.3f", sig_ans(1), sig_ans(2), sig_ans(3));
fclose(f3);

f4 = fopen("H_ans_2.txt", "w");
ans_H_2 = fprintf(f4, "%.3f %.3f %.3f", H_ans_2(1), H_ans_2(2), H_ans_2(3));
fclose(f4);

%% Propagate orbits and attitude with fixed [0.01; -0.01; 0.02] Nm control
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

u_0 = [0.01; -0.01; 0.02]; % Nm, body coords
x_0_att = {I; sigBN_0; omegBN_0; u_0};

out_Attitude = RK4_Attitude(x_0_att, t0, dt, tf, 'torque');

%% Plot simulation outputs
t = out_Attitude(:,1);
sig = out_Attitude(:,2:4);
w = out_Attitude(:,5:7);
u = out_Attitude(:,8:10);

% Orbit
figOrbit = figure;

ax2 = axes();
marsStars = imread("marsStars.jpg");
imshow(marsStars, 'parent', ax2)

ax1 = axes();
title("Orbit Simulation", 'Color', 'w')
hold on
grid on
axis equal
plot3(out_LMO(:,2), out_LMO(:,3), out_LMO(:,4), 'LineWidth', 3)
plot3(out_GMO(:,2), out_GMO(:,3), out_GMO(:,4), 'LineWidth', 3)

mars = surf(marsX, marsY, -marsZ); % Need to flip the sphere for image to map properly
marsSurface = imread("marsSurface.jpg");
set(mars,'FaceColor','texturemap','cdata',marsSurface,'edgecolor','none');

set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w')

view([30 35])

% Attitude
figAttitude = figure;

sgtitle("Nano-satellite Attitude Evolution Over Time")

subplot(3,3,1)
hold on
title("\sigma_1 vs. time")
plot(t, sig(:,1))
xlabel("Time [sec]")
ylabel("\sigma_1")

subplot(3,3,2)
hold on
title("\omega_1 vs. time")
plot(t, w(:,1))
xlabel("Time [sec]")
ylabel("\omega_1")

subplot(3,3,3)
hold on
title("u_1 vs. time")
plot(t, u(:,1))
xlabel("Time [sec]")
ylabel("u_1")

subplot(3,3,4)
hold on
title("\sigma_2 vs. time")
plot(t, sig(:,2))
xlabel("Time [sec]")
ylabel("\sigma_2")

subplot(3,3,5)
hold on
title("\omega_2 vs. time")
plot(t, w(:,1))
xlabel("Time [sec]")
ylabel("\omega_2")

subplot(3,3,6)
hold on
title("u_2 vs. time")
plot(t, u(:,2))
xlabel("Time [sec]")
ylabel("u_2")

subplot(3,3,7)
hold on
title("\sigma_3 vs. time")
plot(t, sig(:,3))
xlabel("Time [sec]")
ylabel("\sigma_3")

subplot(3,3,8)
hold on
title("\omega_3 vs. time")
plot(t, w(:,3))
xlabel("Time [sec]")
ylabel("\omega_3")

subplot(3,3,9)
hold on
title("u_3 vs. time")
plot(t, u(:,3))
xlabel("Time [sec]")
ylabel("u_3")

%% Extract answers to text files
t = 100;

sig_ans_2 = sig(out_Attitude(:,1) == t, :)';

f5 = fopen("sig_ans_2.txt", "w");
ans_sig_2 = fprintf(f5, "%.3f %.3f %.3f", sig_ans_2(1), sig_ans_2(2), sig_ans_2(3));
fclose(f5);


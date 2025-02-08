%% ASEN 5010 Task 9 Main Script
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

% Controller parameters
K = 1/180; % 0.0056, from zeta requirement
P = 1/6; % 0.1667, from time decay requirement
controlParams = [K; P];

% RK4 params
t0 = 0;
dt = 1;
tf = 1000;

%% Propagate orbits and attitude with nadir-pointing control
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

out_Attitude = RK4_Attitude(x_0_att, t0, dt, tf, 'nadir pointing', out_LMO, out_GMO, controlParams);

%% Plot simulation outputs
t = out_Attitude(:,1);
sig = out_Attitude(:,2:4);
w = out_Attitude(:,5:7);
u = out_Attitude(:,8:10);
sigRef = out_Attitude(:,11:13);
wRef = out_Attitude(:,14:16);

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
plot(t, sigRef(:,1), '--')
xlabel("Time [sec]")
ylabel("\sigma_1")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,2)
hold on
grid on
title("\omega_1 vs. time")
plot(t, w(:,1))
plot(t, wRef(:,1), '--')
xlabel("Time [sec]")
ylabel("\omega_1 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,3)
hold on
grid on
title("u_1 vs. time")
plot(t, u(:,1))
xlabel("Time [sec]")
ylabel("u_1 [Nm]")

subplot(3,3,4)
hold on
grid on
title("\sigma_2 vs. time")
plot(t, sig(:,2))
plot(t, sigRef(:,2), '--')
xlabel("Time [sec]")
ylabel("\sigma_2")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,5)
hold on
grid on
title("\omega_2 vs. time")
plot(t, w(:,2))
plot(t, wRef(:,2), '--')
xlabel("Time [sec]")
ylabel("\omega_2 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,6)
hold on
grid on
title("u_2 vs. time")
plot(t, u(:,2))
xlabel("Time [sec]")
ylabel("u_2 [Nm]")

subplot(3,3,7)
hold on
grid on
title("\sigma_3 vs. time")
plot(t, sig(:,3))
plot(t, sigRef(:,3), '--')
xlabel("Time [sec]")
ylabel("\sigma_3")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,8)
hold on
grid on
title("\omega_3 vs. time")
plot(t, w(:,3))
plot(t, wRef(:,3), '--')
xlabel("Time [sec]")
ylabel("\omega_3 [rad/s]")
legend("Current", "Reference", 'location', 'best')

subplot(3,3,9)
hold on
grid on
title("u_3 vs. time")
plot(t, u(:,3))
xlabel("Time [sec]")
ylabel("u_3 [Nm]")

%% Extract answers to text files
time = t;

t = 15;
sig15 = sig(time == t, :);

t = 100;
sig100 = sig(time == t, :);

t = 200;
sig200 = sig(time == t, :);

t = 400;
sig400 = sig(time == t, :);

f1 = fopen("sig15_ans.txt", "w");
ans_sig15 = fprintf(f1, "%.3f %.3f %.3f", sig15(1), sig15(2), sig15(3));
fclose(f1);

f2 = fopen("sig100_ans.txt", "w");
ans_sig100 = fprintf(f2, "%.3f %.3f %.3f", sig100(1), sig100(2), sig100(3));
fclose(f2);

f3 = fopen("sig200_ans.txt", "w");
ans_sig200 = fprintf(f3, "%.3f %.3f %.3f", sig200(1), sig200(2), sig200(3));
fclose(f3);

f4 = fopen("sig400_ans.txt", "w");
ans_sig400 = fprintf(f4, "%.3f %.3f %.3f", sig400(1), sig400(2), sig400(3));
fclose(f4);




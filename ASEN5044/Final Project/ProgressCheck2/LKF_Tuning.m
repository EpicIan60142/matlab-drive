%% ASEN 5044 Project 2 Progress Check LKF Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
load("orbitdeterm_finalproj_KFdata.mat");

r0 = 6678; % km
muEarth = 398600; % km^3/s^2
n = sqrt(muEarth/(r0^3)); % rad/s

dt = 10; % sec

t_l = tvec;

X0 = r0; % km
n = sqrt(muEarth/(r0^3)); % rad/s
YDot0 = r0*n; % km/s

xNom0 = [X0; 0; 0; YDot0];
xPerturb0 = [0; 0.075; 0; -0.021];
x0 = xNom0 + 0*xPerturb0;

P0 = 10*eye(size(xPerturb0,1));

Q_LKF = eye(4); %Qtrue;
% Q_LKF = 1000*diag([1000 100 1000 100]);
R_LKF = Rtrue;

options = odeset('RelTol',10^(-12), 'AbsTol',10^(-12));
[~, xNom_l] = ode45(@(t, x) nonlinearOrbitalModel(x, [0;0], [0;0], muEarth, t), t_l, xNom0, options);
xNom_l = xNom_l';

filterParams = struct('Q', Q_LKF, 'R', R_LKF);
filterInit = struct('dx0', xPerturb0, 'P0', P0, 'xNom', xNom_l, 'dt', dt);

%% Run LKF
LKF_out = LKF(ydata, filterParams, filterInit, [0;0]);
t_KF = LKF_out.t_KF;
x_KF = LKF_out.x_KF;

% Plot trajectory
figure
hold on; grid on;
plot(xNom_l(1,:), xNom_l(3,:), 'b-')
plot(x_KF(1,:), x_KF(3,:), 'r--');

%% Run NEES/NIS
filterFunc = @LKF;
N = 1;
alpha = 0.05;
timesteps = 500; %tvec(end)/dt;

[epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(filterFunc, N, alpha, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0);

%% Plot NEES/NIS
figure
hold on; grid on;
title("NEES Statistic")
plot(epsilon_x_bar, 'ro');
yline(NEESbounds, 'r--');

figure
hold on; grid on;
title("NIS Statistic")
plot(epsilon_y_bar, 'bo');
plot(NISbounds(1,:), 'b--')
plot(NISbounds(2,:), 'b--')

%% Plot errors
t_KF = LKF_out.t_KF;
x_KF = LKF_out.x_KF;
sig_KF = LKF_out.sig_KF;

stateError = xNom_l - x_KF;

figure
ax(1) = subplot(2,2,1);
    hold on; grid on;
    title("X state error")
    plot(t_KF, stateError(1,:))
    plot(t_KF, 2*sig_KF(1,:), 'r--');
    plot(t_KF, -2*sig_KF(1,:), 'r--')
    xlabel("Time [sec]"); ylabel("X [km]")
ax(2) = subplot(2,2,2);
    hold on; grid on;
    title("Xdot state error")
    plot(t_KF, stateError(2,:))
    plot(t_KF, 2*sig_KF(2,:), 'r--');
    plot(t_KF, -2*sig_KF(2,:), 'r--')
ax(3) = subplot(2,2,3);
    hold on; grid on;
    title("Y state error")
    plot(t_KF, stateError(3,:))
    plot(t_KF, 2*sig_KF(3,:), 'r--');
    plot(t_KF, -2*sig_KF(3,:), 'r--')
ax(4) = subplot(2,2,4);
    hold on; grid on;
    title("Ydot state error")
    plot(t_KF, stateError(4,:))
    plot(t_KF, 2*sig_KF(4,:), 'r--');
    plot(t_KF, -2*sig_KF(4,:), 'r--')
linkaxes(ax, 'x');



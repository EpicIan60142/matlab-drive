clc 
clear
close all

options = odeset('RelTol',1e-12, "AbsTol",1e-12);

%time vars
totalTime = 14000;
dT = 10;
thist = 0:dT:14000;

%constants
mu = 398600;
rE = 6378;
omegaE = 2*pi/86400;
r0 = 6678;
n = sqrt(mu/(r0^3));
T = 2*pi/n;

%initial states
perturbation = [0;0.075;0;-0.021];
x0 = [6678; 0;0; 6678*sqrt(mu/(6678^3))] + 1*perturbation;

%nominal solution
Xnom = r0*cos(sqrt(mu/r0^3).*thist);
Xdotnom = -sqrt(mu/r0^3)*r0*sin(sqrt(mu/r0^3).*thist);
Ynom = r0*sin(sqrt(mu/r0^3).*thist);
Ydotnom = r0*sqrt(mu/r0^3)*cos(sqrt(mu/r0^3).*thist);

nomSolution = [Xnom; Xdotnom; Ynom; Ydotnom];


% Assign input structures and matrices
load("orbitdeterm_finalproj_KFdata.mat")
P0 = 10.*eye(4);

Gamma = [0 0; 1 0; 0 0; 0 1];


%constants
consts(1) = rE;
consts(2) = omegaE;
consts(3) = mu;

%LKF
Q_LKF = 1e-7*[
                5       1e-3    1e-3      0
                1e-3    1e-1    0        1e-3
                1e-3    0       5        1e-3
                0       1e-3    1e-3     1e-1
             ];
filterParams = struct('R', Rtrue, 'Q', Q_LKF); %1e-7*diag([0.005 10 0.005 10]));
filterInit = struct('dx0', perturbation, 'P0', P0, 'xNom', nomSolution, "dt", dT);

%EKF
if false
    filterParams.Gamma = Gamma;
    filterParams.Q = eye(4);
    filterParams.R = Rtrue;
    filterParams.dT = dT;
    filterParams.tvec = tvec;
    filterParams.consts = consts;
end
%initialize filter
perturb_x0 = [0;0.075;0;-0.021];
xhat0 = x0;
P0 = 1e6*eye(4);

P0 = diag([5 1 5 1]);

filterInit.xhat0 = xhat0;
filterInit.P0 = P0;

[epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(@LKF, 100, 0.05, filterParams, Qtrue, Rtrue, filterInit, 0*50 + 1*floor(T/dT), x0);
% [epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds] = chi2tests(@EKF_Func, 50, 0.05, filterParams, Qtrue, Rtrue, filterInit, 400, x0);

% EKFdata = EKF_Func(ydata, filterParams, filterInit);

figure()
hold on; grid on;
title("NEES Statistic")
plot(epsilon_x_bar, "ro")
yline(NEESbounds(1), "k--")
yline(NEESbounds(2), "k--")

figure()
hold on; grid on;
title("NIS Statistic")
plot(epsilon_y_bar, "bo")
plot(NISbounds(1,:), "k--")
plot(NISbounds(2,:), "k--")
ylim([0,10])

t_KF = LKF_out.t_KF;
x_KF = LKF_out.x_KF;
sig_KF = LKF_out.sig_KF;

stateError = nomSolution(:,1:size(x_KF,2)) - x_KF;

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

% figure
% hold on;
% plot(nomSolution(1,1:size(x_KF,2)), nomSolution(3,1:size(x_KF,2)), 'b-')
% plot(x_KF(1,:), x_KF(3,:), 'r--');





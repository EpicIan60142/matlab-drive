%% ASEN 5044 Project 2 Progress Check LKF Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
load("orbitdeterm_finalproj_KFdata.mat");

    % Physical constants
muEarth = 398600; % km^3/s^2
rEarth = 6378; % km
omegaEarth = (2*pi/86400); % rad/s

    % Satellite
r0 = 6678; % km
X0 = r0; % km
n = sqrt(muEarth/(r0^3)); % rad/s
YDot0 = r0*n; % km/s

xNom0 = [X0; 0; 0; YDot0];
xPerturb0 = [0; 0.075; 0; -0.021];
x0 = xNom0 + xPerturb0;

P0 = 999*eye(size(xNom0,1));

    % Filter
dT = 10; % sec

%% Test LKF
    % Define Q and R for the LKF
Q_LKF = 0.01*Qtrue;
R_LKF = 1000*Rtrue;

    % Calculate nominal trajectory and station states
t_l = tvec;

xNom_l = zeros(size(xNom0,1), length(t_l));
for k = 1:length(t_l)
    xNom_l(:,k) = [
                        r0*cos(n*t_l(k));
                        -r0*n*sin(n*t_l(k));
                        r0*sin(n*t_l(k));
                        r0*n*cos(n*t_l(k))
                  ];
end

    % Assign input structures and matrices
filterParams = struct('Q', Q_LKF, 'R', R_LKF);
filterInit = struct('dx0', xPerturb0, 'P0', P0, 'xNom', xNom_l, 'dt', dT);
debug = [0; 0];

    % Run LKF
LKF_out = LKF(ydata, filterParams, filterInit, debug);

%% Plot results
t_KF = LKF_out.t_KF;
x_KF = LKF_out.x_KF;
sig_KF = LKF_out.sig_KF;

stateError = xNom_l - x_KF;

figure
hold on; grid on; axis equal
title("TMT vs. KF Trajectory")
plot(xNom_l(1,:), xNom_l(3,:),'b-')
plot(x_KF(1,:), x_KF(3,:),'r--')
xlabel("X [km]"); ylabel("Y [km]")
legend("TMT", "LKF")

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


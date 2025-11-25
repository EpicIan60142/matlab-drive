%% ASEN 6014 Coursera Relative Motion Control Quiz 15 Main
% By: Ian Faber, 11/16/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Planetary constants structure
pConst.mu = 398600.4415; % km^3/s^2

    % Chief initial orbit elements
oe_c.mu = pConst.mu;
oe_c.a = 7500; % km
oe_c.e = 0;
oe_c.i = deg2rad(0); % deg -> rad
oe_c.RAAN = deg2rad(0); % deg -> rad
oe_c.argPeri = deg2rad(0); % deg -> rad
oe_c.f = deg2rad(0); % deg -> rad

    % Desired CWH coefficients
coeff.A0 = 2; % km
coeff.B0 = 0.5; % km
coeff.xOff = 0; % km
coeff.yOff = 10; % km
coeff.alpha = 0;
coeff.beta = 0;

    % Control gains
kConst.K = 1e-5;
kConst.P = 4e-3;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Run controller
    % Setup initial state vector
cart_c = convClassicOE2Cart(oe_c);
X0_c = [cart_c.rVec; cart_c.vVec];

X0_d = [1; -2; 2; 0; 0; 0];

X0 = [X0_c; X0_d];

    % Define tspan to output at specified times 0s, 1000s, and 4848s
tspan = 0:10:1*4500;

    % Run controller
[t, X] = ode45(@(t,X)LVLHFeedbackControlCoeffs(t,X,coeff,kConst,pConst), tspan, X0, opt);

[~, u, X_r] = cellfun(@(t,X)LVLHFeedbackControlCoeffs(t,X.', coeff, kConst, pConst), num2cell(t), num2cell(X,2), 'uni', 0);
u = cell2mat(cellfun(@(x)x',u,'uni',0));
X_r = cell2mat(cellfun(@(x)x', X_r,'uni',0));

%% Plot and report undisturbed results
    % Extract states
X_c = X(:,1:6);
X_d_Hill = X(:,7:12);

    % Plot
markerSize = 30;
figure; axis equal
hold on; grid on;
title("LVLH Relative Motion Control")
scatter3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), markerSize, 'g', 'filled');
plot3(X_d_Hill(:,1), X_d_Hill(:,2), X_d_Hill(:,3), 'b-')
plot3(X_r(:,1), X_r(:,2), X_r(:,3), 'k--')
scatter3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), markerSize, 'r', 'filled')
scatter3(0,0,0,markerSize,'k','*')
xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
view([30,35]); xlim([-2.5, 2.5]); ylim([-5, 15]); zlim([-3, 2.5]);
legend("Start", "Trajectory", "Desired", "End", "Chief")

figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Controls")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,1), 'b-');
    xlabel("Time [sec]"); ylabel("u_1 [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,2), 'b-');
    xlabel("Time [sec]"); ylabel("u_2 [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,3), 'b-');
    xlabel("Time [sec]"); ylabel("u_3 [m/s^2]");
linkaxes(ax, 'x');

figure; tl = tiledlayout(3,2); ax = [];
title(tl, "Differences")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,1) - X_r(:,1), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta x [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,2) - X_r(:,2), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta y [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,3) - X_r(:,3), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta z [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,4) - X_r(:,4), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta xDot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,5) - X_r(:,5), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta yDot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X_d_Hill(:,6) - X_r(:,6), 'k--')
    xlabel("Time [sec]"); ylabel("\Delta zDot [m/s]");
linkaxes(ax, 'x');

    % Report
fprintf("rho0 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==0, 1:3));
fprintf("rho1200 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==1200, 1:3));
fprintf("rho4500 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==4500, 1:3));




%% ASEN 6014 HW 1 Coursera Quiz 5 Question 4 Main Script
% By: Ian Faber, 08/26/2025

%% Housekeeping
clc; clear; close all;

%% Setup
    % Satellite 1
r1_init = [-6685.20926; 601.51244; 3346.06634]; % km
v1_init = [-1.74294; -6.70242; -2.27739]; % km/s
X1_init = [r1_init; v1_init];

    % Satellite 2
r2_init = [-6685.21657; 592.52839; 3345.6716]; % km
v2_init = [-1.74283; -6.70475; -2.27334]; % km/s
X2_init = [r2_init; v2_init];

X0 = [X1_init; X2_init];

t_end = 4848; % 4848 seconds of propagation time

tspan = 0:t_end;

J2_enable = true;

pConst.mu = 398600.4415; % km^3/s^2
pConst.Ri = 6378; % km
pConst.J2 = 1.08263e-3; % n.d.

opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% Propagate states
[t, X] = ode45(@(t,X)multiSatOrbitEOM(t, X, pConst, J2_enable), tspan, X0, opt);

%% Pull out states at the end of propagation
format long

r1_f = X(end, 1:3)'
v1_f = X(end, 4:6)'

r2_f = X(end, 7:9)'
v2_f = X(end, 10:12)'

    % Coursera answer for convenience
r1 = [1752.88687, 6891.31557, 2375.24597];
v1 = [-6.51517, 0.54618, 3.22361];
r2 = [1760.55329, 6882.87250, 2362.55702];
v2 = [-6.51985, 0.55669, 3.24430];

%% Plot states
figure;
hold on; grid on; axis equal;
title("Integrated Spacecraft Motion")
plot3(X(:,1), X(:,2), X(:,3), 'b-');
plot3(X(:,7), X(:,8), X(:,9), 'r-');
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
view([30 35]); legend("Spacecraft 1", "Spacecraft 2")

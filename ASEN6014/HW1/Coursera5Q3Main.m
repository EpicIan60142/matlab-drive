%% ASEN 6014 HW 1 Coursera Quiz 5 Question 3 Main Script
% By: Ian Faber, 08/26/2025

%% Housekeeping
clc; clear; close all;

%% Setup
r0 = [2466.69; 5941.54; 3282.71]; % km
v0 = [-6.80822; 1.04998; 3.61939]; % km/s
X0 = [r0; v0];

t_end = 60*60; % 60 minutes of propagation time

tspan = 0:t_end;

J2_enable = false;

pConst.mu = 398600.4415; % km^3/s^2
pConst.Ri = 6378; % km
pConst.J2 = 1.08264e-3; % n.d.

opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% Propagate states
[t, X] = ode45(@(t,X)orbitEOM(t, X, pConst, J2_enable), tspan, X0, opt);

%% Pull out states at the end of propagation
rf = X(end, 1:3)'
vf = X(end, 4:6)'

    % Coursera answer for convenience
rVec_N = [-3437.2, -7152.8, -3755.7];
vVec_N = [5.5743, -0.9220, -3.0087];

%% Plot integrated states
figure;
hold on; grid on; axis equal;
title("Integrated Spacecraft Motion")
plot3(X(:,1), X(:,2), X(:,3), 'b-');
xlabel("X [km]"); ylabel("Y [km]"); zlabel("Z [km]");
view([30 35])





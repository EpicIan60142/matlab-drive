%% ASEN 6014 Formation Flying Coursera Quiz 6 Main Script
% By: Ian Faber, 09/16/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Earth mu
const.mu = 398600.4415; % km^3/s^2

    % Chief initial state
oe.mu = const.mu;
oe.a = 6800;
oe.e = 0;
oe.i = deg2rad(5);
oe.RAAN = deg2rad(10);
oe.argPeri = deg2rad(15);
oe.f = deg2rad(20);
cart = convClassicOE2Cart(oe);
X_c_N = [cart.rVec; cart.vVec];

    % Deputy initial state
r_d_H = [1.299038; -1.0000; 0.3213938]; % km
v_d_H = [-0.000844437; -0.00292521; -0.000431250]; % km/s
X_d_H = [r_d_H; v_d_H];

    % Integration time and options
tspan = 0:10:1300; % sec
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Task 1: Integrate both chief and deputy in inertial frame
    % Convert deputy to inertial
X0_N = [X_c_N; convDeputyH2N(X_c_N, X_d_H, const)];

    % Integrate in the inertial frame
[t_N, X_N] = ode45(@(t,X)multiSatOrbitEOM(t,X,const,false), tspan, X0_N, opt);

    % Pull out final deputy states
fprintf("\nInertial - Deputy final inertial state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_N(end, 7:12))

%% Task 2: Integrate relative equations of motion
    % Populate initial states
X0 = [X_c_N; X_d_H];

    % Integrate full nonlinear relative EOM
[t_rel, X_rel] = ode45(@(t,X)fullCWHEOM(t,X,const), tspan, X0, opt);

    % Pull out final deputy states
fprintf("\nRelative - Deputy final inertial state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", convDeputyH2N(X_rel(end,1:6)', X_rel(end, 7:12)', const))
fprintf("\nDeputy final Hill state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_rel(end, 7:12))

%% Task 3: Investigate how much a deputy shifts compared to what CWH would predict
    % Set semi-major axis to 7500 km
oe.a = 7500;
oe.e = 0;
oe.i = 0;
oe.RAAN = 0;
oe.argPeri = 0;
oe.f = 0;
cart = convClassicOE2Cart(oe);
X_c_N = [cart.rVec; cart.vVec];

    % Initialize deputy with along-track offset of 200 km only
r_d_H = [0; 200; 0];
v_d_H = [0; 0; 0];
X_d_H = [r_d_H; v_d_H];
X_d_N = convDeputyH2N(X_c_N, X_d_H, const);

    % Initial states
X0_N = [X_c_N; X_d_N];
X0 = [X_c_N; X_d_H];

    % Integrate longer
tspan = 0:10:2000; % sec

fprintf("\nDeputy initial relative state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_d_H)

    % Integrate full nonlinear EOM
[t_N, X_N] = ode45(@(t,X)multiSatOrbitEOM(t,X,const,false), tspan, X0_N, opt);
X_rel = convDeputyN2H(X_N(end,1:6)', X_N(end,7:12)', const);

[t_rel, X_rel] = ode45(@(t,X)fullNLGenRelEOM(t,X,const), tspan, X0, opt);
X_rel = X_rel(end, 7:12)';

fprintf("\nDeputy final relative state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_rel)

fprintf("\nFinal separation distance: %.6f km\n", norm(X_rel(1:3)))

fprintf("\nChange in separation distance: %.4f km\n", norm(X_rel(1:3)) - norm(r_d_H))

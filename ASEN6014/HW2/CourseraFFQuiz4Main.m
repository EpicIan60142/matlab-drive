%% ASEN 6014 Formation Flying Coursera Quiz 4 Main Script
% By: Ian Faber, 09/16/2025

%% Housekeeping
clc; clear; close all;

%% Setup
    % Earth mu
const.mu = 398600.4415; % km^3/s^2

sc = zeros(2,1);

    % Chief initial state
r_c_N = [-6685.20926; 601.51244; 3346.06634]; % km
v_c_N = [-1.74294; -6.70242; -2.27739]; % km/s
X_c_N = [r_c_N; v_c_N];

    % Deputy initial state
r_d_H = [-81.22301; 248.14201; 94.95904]; % km
v_d_H = [0.47884; 0.14857; 0.13577]; % km/s
X_d_H = [r_d_H; v_d_H];

    % Integration time and options
tspan = 0:10:2000;
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);

%% Task 1: Integrate both chief and deputy in inertial frame
    % Convert deputy to inertial
X0_N = [X_c_N; convDeputyH2N(X_c_N, X_d_H, const)];
fprintf("\nDeputy initial inertial state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X0_N(7:12))

    % Integrate in the inertial frame
[t_N, X_N] = ode45(@(t,X)multiSatOrbitEOM(t,X,const,sc,false), tspan, X0_N, opt);

    % Pull out final deputy states
fprintf("\nInertial - Deputy final inertial state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_N(end, 7:12))

%% Task 2: Integrate relative equations of motion
    % Populate initial states
X0 = [X_c_N; X_d_H];

    % Integrate full nonlinear relative EOM
[t_rel, X_rel] = ode45(@(t,X)fullNLGenRelEOM(t,X,const), tspan, X0, opt);

    % Pull out final deputy states
fprintf("\nRelative - Deputy final inertial state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", convDeputyH2N(X_rel(end,1:6)', X_rel(end, 7:12)', const))
fprintf("\nDeputy final Hill state: [%.6f, %.6f, %.6f] km, [%.6f, %.6f, %.6f] km/s\n", X_rel(end, 7:12))



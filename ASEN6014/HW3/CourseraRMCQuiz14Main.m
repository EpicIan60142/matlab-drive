%% ASEN 6014 HW 3 Coursera Relative Motion Control Quiz 14 Main
% By: Ian Faber, 11/16/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
format long

    % Planetary constant
mu = 398600.4415; % km^3/s^2

    % Deputy elements
a_d = 8500; % km
e_d = 0.1;

    % Desired differences
deltaE = 0.002;
deltaA = 0;

    % Intermediate values
eta = sqrt(1-e_d^2);
n = sqrt(mu/(a_d^3));
h = sqrt(mu*a_d*eta^2);
r_p = a_d*(1-e_d);
r_a = a_d*(1+e_d);

%% Calculate delta V's
dVThetaP = ((n*a_d*eta)/4)*(deltaA/a_d + deltaE/(1+e_d))

dVThetaA = ((n*a_d*eta)/4)*(deltaA/a_d - deltaE/(1-e_d))

%% Execute burn at periapsis first
fprintf("Correction at periapsis first:\n")
    % Original and new periapsis velocities
v_p0 = h/r_p;
v_p1 = v_p0 + dVThetaP;

    % New semi major axis
a_1 = -mu/(v_p1^2 - (2*mu/r_p));

    % Calculate anomalies
M_0 = convE2M(convf2E(0, e_d), e_d);
M_f = convE2M(convf2E(pi, e_d), e_d);

    % Calculate change in anomaly
da = a_1 - a_d;

dM_rad = -1.5*(M_f-M_0)*((a_1-a_d)/a_d)
dM_deg = rad2deg(dM_rad)

%% Execute burn at apoapsis first
fprintf("Correction at apoapsis first:\n")
    % Original and new apoapsis velocities
v_a0 = h/r_a;
v_a1 = v_a0 + dVThetaA;

    % New semi major axis
a_1 = -mu/(v_a1^2 - (2*mu/r_a));

    % Calculate anomalies
M_0 = convE2M(convf2E(0, e_d), e_d);
M_f = convE2M(convf2E(pi, e_d), e_d);

    % Calculate change in anomaly
da = a_1 - a_d;

dM_rad = -1.5*(M_f-M_0)*((a_1-a_d)/a_d)
dM_deg = rad2deg(dM_rad)

format default


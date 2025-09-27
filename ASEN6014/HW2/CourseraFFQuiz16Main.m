%% ASEN 6014 Coursera Formation Flying Quiz 16 Main Script
% By: Ian Faber, 09/26/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Chief spacecraft
oe_c.a = 10000; % km
oe_c.e = 0.2;
oe_c.i = deg2rad(37);
oe_c.RAAN = deg2rad(40);
oe_c.argPeri = deg2rad(65);
oe_c.f = deg2rad(10);

    % Deputy element differences
doe.da = 0; % km
doe.de = 0.0001;
doe.di = 0.001; % rad
doe.dRAAN = 0.001; % rad
doe.dargPeri = 0.001; % rad
doe.dM0 = -0.001; % rad

%% Calculate solutions at specific true anomalies
f1 = oe_c.f;
rho1 = calcLinearGenEllSol(oe_c, doe, f1);

fprintf("rho_0_H = [%.6e, %.6e, %.6e]\n", rho1);

f2 = oe_c.f + deg2rad(60);
rho2 = calcLinearGenEllSol(oe_c, doe, f2);

fprintf("rho_t_H = [%.6e, %.6e, %.6e]\n", rho2);


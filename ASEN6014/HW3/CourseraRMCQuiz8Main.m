%% ASEN 6014 Coursera Relative Motion Control Quiz 8 Main
% By: Ian Faber, 11/07/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Planetary constants structure
pConst.mu = 398600.4415; % km^3/s^2

    % Chief initial orbit elements
oe_c.mu = pConst.mu;
oe_c.a = 8500; % km
oe_c.e = 0.2;
oe_c.i = deg2rad(53); % deg -> rad
oe_c.RAAN = deg2rad(55); % deg -> rad
oe_c.argPeri = deg2rad(40); % deg -> rad
oe_c.f = deg2rad(0); % deg -> rad

    % Initial deputy orbit element differences
doe.de = 1e-3;
doe.di = deg2rad(0.1);
doe.dRAAN = deg2rad(0.1);

    % Initial deputy orbit elements
oe_d = oe_c;
oe_d.e = oe_c.e + doe.de;
oe_d.i = oe_c.i + doe.di;
oe_d.RAAN = oe_c.RAAN + doe.dRAAN;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);



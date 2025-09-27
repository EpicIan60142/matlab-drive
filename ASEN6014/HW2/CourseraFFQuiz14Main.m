%% ASEN 6014 Coursera Formation Flying Quiz 14 Main Script
% By: Ian Faber, 09/26/2025

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
    % Earth gravitational parameter
mu = 398600.4415;

    % Chief orbital elements
oe_c.mu = mu;
oe_c.a = 7000; % km
oe_c.e = 0.01;
oe_c.i = deg2rad(78);
oe_c.RAAN = deg2rad(120);
oe_c.argPeri = deg2rad(33);
oe_c.f = deg2rad(45);

    % Deputy orbital elements
oe_d = oe_c;
oe_d.f = oe_c.f + deg2rad(1);

%% Calculate cartesian center of mass
    % Calculate cartesian vectors
cart_c = convClassicOE2Cart(oe_c);
r_c = cart_c.rVec;
v_c = cart_c.vVec;

cart_d = convClassicOE2Cart(oe_d);
r_d = cart_d.rVec;
v_d = cart_d.vVec;

    % Combine vectors for the system
r_system = [r_c, r_d];
v_system = [v_c, v_d];

    % Calculate center of mass, assuming identical masses for spacecraft
r_center = zeros(size(r_c));
for k = 1:size(r_system,2)
    r_center = r_center + r_system(:,k);
end
r_center = r_center/size(r_system,2);

fprintf("rCCM_N = [%.6e,%.6e,%.6e]\n", r_center)

v_center = zeros(size(v_c));
for k = 1:size(v_system,2)
    v_center = v_center + v_system(:,k);
end
v_center = v_center/size(v_system,2);

%% Calculate orbit element center of mass
    % Combine all orbit elements into one matrix
oe_system = [[oe_c.a; oe_c.e; oe_c.i; oe_c.RAAN; oe_c.argPeri; oe_c.f],...
            [oe_d.a; oe_d.e; oe_d.i; oe_d.RAAN; oe_d.argPeri; oe_d.f]];

    % Calculate center of mass, assuming identical spacecraft masses
oe_center = zeros(size(oe_system,1));
for k = 1:size(oe_system,2)
    oe_center = oe_center + oe_system(:,k);
end
oe_center = oe_center/size(oe_system,2);

    % Convert result to cartesian
oe_cent.mu = mu;
oe_cent.a = oe_center(1);
oe_cent.e = oe_center(2);
oe_cent.i = oe_center(3);
oe_cent.RAAN = oe_center(4);
oe_cent.argPeri = oe_center(5);
oe_cent.f = oe_center(6);

cart_cent = convClassicOE2Cart(oe_cent);
r_cent = cart_cent.rVec;

fprintf("rOECM_N = [%.6e,%.6e,%.6e]\n", r_cent)

%% Convert the cartesian center of mass classical orbit elements and determine the difference in orbit periods
cart_center.mu = mu;
cart_center.rVec = r_center;
cart_center.vVec = v_center;
oe_center = convCart2ClassicOE(cart_center);

a_center = oe_center.a;

P_center = 2*pi*sqrt((a_center^3)/mu);
P_OEcenter = 2*pi*sqrt((oe_cent.a^3)/mu);

Pdiff = abs(P_OEcenter - P_center)

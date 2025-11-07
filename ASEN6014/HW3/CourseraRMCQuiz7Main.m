%% ASEN 6014 Coursera Relative Motion Control Quiz 7 Main
% By: Ian Faber, 11/07/2025

%% Housekeeping
clc; clear; close all;

addpath("..\utilities\")

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
doe.e = 1e-3;
doe.i = deg2rad(0.1);
doe.RAAN = deg2rad(0.1);

    % Initial deputy orbit elements
oe_d = oe_c;
oe_d.e = oe_c.e + doe.e;
oe_d.i = oe_c.i + doe.i;
oe_d.RAAN = oe_c.RAAN + doe.RAAN;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Integrate chief and deputy forward in time
    % Set up initial state vector
cart_c = convClassicOE2Cart(oe_c);
X0_c = [cart_c.rVec; cart_c.vVec];

cart_d = convClassicOE2Cart(oe_d);
X0_d = [cart_d.rVec; cart_d.vVec];
X0_d_H = convDeputyN2H(X0_c, X0_d, pConst);

X0 = [X0_c; X0_d_H];

    % Define tspan to output at specified times 0s, 1000s, and 4000s
tspan = 0:10:8000;

    % Integrate relative motion
[t, X] = ode45(@(t,X)fullNLGenRelEOM(t,X,pConst), tspan, X0, opt);

%% Plot and report reports
    % Extract results
X_d_Hill = X(:,7:12);

    % Plot
figure; axis equal
hold on; grid on;
title("OE-based Reference Trajectory")
plot3(X_d_Hill(:,1), X_d_Hill(:,2), X_d_Hill(:,3), 'b-')
xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
view([30,35]);

    % Report
fprintf("rho0 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==0, 1:3));
fprintf("rho1000 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==1000, 1:3));
fprintf("rho4000 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==4000, 1:3));






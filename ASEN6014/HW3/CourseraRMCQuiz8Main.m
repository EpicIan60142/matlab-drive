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
doe.da = 1; % km
doe.de = 1e-3;
doe.di = deg2rad(0.1);
doe.dRAAN = deg2rad(0.1);

    % Initial deputy orbit elements
oe_d = oe_c;
oe_d.a = oe_c.a + doe.da;
oe_d.e = oe_c.e + doe.de;
oe_d.i = oe_c.i + doe.di;
oe_d.RAAN = oe_c.RAAN + doe.dRAAN;

    % Desired orbit element differences
doe_r.da = 0;
doe_r.de = 5e-4;

    % Controller parameters
        % P11
kConst.P11Off = 1e-12;
kConst.P11Amp = 1e-3;
kConst.P11Exp = 10;

        % P22
kConst.P22Off = 1e-12;
kConst.P22Amp = 1e-3;
kConst.P22Exp = 10;

    % ode45 settings
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Run controller
    % Setup initial state vector
cart_c = convClassicOE2Cart(oe_c);
X0_c = [cart_c.rVec; cart_c.vVec];

cart_d = convClassicOE2Cart(oe_d);
X0_d = [cart_d.rVec; cart_d.vVec];

X0 = [X0_c; X0_d];

    % Define tspan to output at specified times 0s, 1000s, and 4848s
tspan = 0:10:28000;

    % Run controller
[t, X] = ode45(@(t,X)OEFeedbackControlElems(t,X,doe_r,kConst,pConst), tspan, X0, opt);

%% Plot and report results
    % Extract states
X_c = X(:,1:6);
X_d = X(:,7:12);

    % Convert to Hill Frame
X_d_Hill = zeros(size(X_d));
for k = 1:length(t)
    X_d_Hill(k,:) = convDeputyN2H(X_c(k,:)', X_d(k,:)', pConst)';
end

    % Plot
markerSize = 30;
figure; axis equal
hold on; grid on;
title("OE Difference Relative Motion Control")
scatter3(X_d_Hill(1,1), X_d_Hill(1,2), X_d_Hill(1,3), markerSize, 'g', 'filled');
plot3(X_d_Hill(:,1), X_d_Hill(:,2), X_d_Hill(:,3), 'b-')
scatter3(X_d_Hill(end,1), X_d_Hill(end,2), X_d_Hill(end,3), markerSize, 'r', 'filled')
% scatter3(0,0,0,markerSize,'k','*')
xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]")
view([30,35]);
legend("Start", "Trajectory", "End")%, "Chief")

    % Report
fprintf("rho0 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==0, 1:3));
fprintf("rho1600 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==1600, 1:3));
fprintf("rho7800 = [%.3e, %.3e, %.3e]\n", X_d_Hill(t==7800, 1:3));


%% ASEN 6020 HW 1 Problem 6 Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
p = 10000; % km
mu = 398600.4415; % km^3/s^2 

dV_a = @(e) 2*sqrt(mu/p)*e;
dV_b = @(e) 2*sqrt(mu/p)*(sqrt(1-e) - (1-e));
dV_c = @(e) 2*sqrt(mu/p)*(1 + e - sqrt(1+e));

%% Compare costs over the valid range of eccentricity for this problem
e = 0:0.001:1;

figure
hold on; grid on;
titleText = sprintf("Argument of Periapsis Rotation Maneuver Cost Comparison,\ngiven p = %.1f km and \\mu = %.4f km^3/s^2", p, mu);
title(titleText)
plot(e, dV_a(e), 'b-');
plot(e, dV_b(e), 'r-');
plot(e, dV_c(e), 'k-');
xlabel("Eccentricity"); ylabel("Cost (\DeltaV) [km/s]");
legend("One-Impulse", "Two-Impulse, Apoapsis", "Two-Impulse, Periapsis", 'Location','best')

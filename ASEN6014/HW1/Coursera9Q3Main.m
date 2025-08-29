%% ASEN 6014 HW 1 Coursera Quiz 9 Question 3 Main Script
% By: Ian Faber, 08/27/2025

%% Housekeeping
clc; clear; close all;

addpath('../utilities/')

%% Setup
mu = 398600.4415; % km^3/s^2
a = 7500; % km
e = 0.05;
f0 = deg2rad(25); % rad

t = 3600; % seconds

%% Solve for M 1 hour later
E0 = convf2E(f0,e);
M0 = E0 - e*sin(E0);

n = sqrt(mu/(a^3));

M = M0 + n*(t-0);

%% Convert to true anomaly
E = convM2E(M,e);

f = convE2f(E,e); % rad

if f < 0
    f = f + 2*pi
else
    f
end
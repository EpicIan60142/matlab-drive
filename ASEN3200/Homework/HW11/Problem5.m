%% ASEN 3200 HW O-5 Problem 5 Script
%   By: Ian Faber, 11/29/2022

%% Housekeeping
clc; clear; close all

%% Constants and givens
mu = 398600; % km^3/s^2, Earth
rA = 300 + 6378; % km
T = 2*pi*sqrt(rA^3/mu); % sec

n = (2*pi)/T; % rad/s

r0 = [50; 0; 0]; % m
v0 = [0; -100*n; 0]; % m/s

size = 35; % Marker size for 3D plot

%% Relative Motion matrices
Phi_rr = @(t,n) [
                    4 - 3*cos(n*t),         0,      0;
                    6*(sin(n*t) - n*t),     1,      0;
                    0,                      0,      cos(n*t)
                ];

Phi_rv = @(t,n) [
                    (1/n)*sin(n*t),         (2/n)*(1-cos(n*t)),     0;
                    (2/n)*(cos(n*t) - 1),   (4/n)*sin(n*t) - 3*t,   0;
                    0,                      0,                      (1/n)*sin(n*t)
                ];

Phi_vr = @(t,n) [
                    3*n*sin(n*t),           0,      0;
                    6*(n*cos(n*t)-n),       0,      0;
                    0,                      0,      -n*sin(n*t)
                ];

Phi_vv = @(t,n) [
                    cos(n*t),       2*sin(n*t),         0;
                    -2*sin(n*t),    4*cos(n*t) - 3,     0;
                    0,              0,                  cos(n*t)
                ];


%% Problem

v0plus = (Phi_rv(T/4,n)^-1)*([0; -200; 0] - Phi_rr(T/4, n)*[50; 0; 0]); % m/s
deltaV1Vec = v0plus - v0
deltaV1 = norm(deltaV1Vec) % m/s

v1minus = Phi_vr(T/4,n)*r0 + Phi_vv(T/4,n)*v0plus
deltaV2Vec = [0; 0; 0] - v1minus
deltaV2 = norm(deltaV2Vec)





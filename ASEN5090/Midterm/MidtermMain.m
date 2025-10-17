%% ASEN 5090 Midterm Main Script
% By: Ian Faber, 10/16/2026

%% Housekeeping
clc; clear; close all;

addpath("../utilities/")

%% Setup
mu = 398600.5; % km^3/s^2
c = 3e5; % km/s

%% Problem 1
    % Setup
r_1_ECEF = [1193.698; -4454.942; 26156.494]; % km
r_2_ECEF = [832.273; -4720.052; 4792.867]; % km
r_3_ECEF = [-3994.192; 22652.188; 13280]; % km
r_4_ECEF = [21189.069; -36700.543; 0]; % km
r_5_ECEF = [7867.022; -21614.466; 13280]; % km

r_A_LLA = [26, -80, 0];

    % Part a
v1 = sqrt(mu/norm(r_1_ECEF))
v2 = sqrt(mu/norm(r_2_ECEF))
v3 = sqrt(mu/norm(r_3_ECEF))
v4 = sqrt(mu/norm(r_4_ECEF))
v5 = sqrt(mu/norm(r_5_ECEF))

    % Part b
format long
r_A_ECEF = lla2ecef(r_A_LLA, 'WGS84')' % m
format default

r_A_ECEF = r_A_ECEF/1000; % m -> km

    % Part c
satECEF = [r_1_ECEF, r_2_ECEF, r_3_ECEF, r_4_ECEF, r_5_ECEF];
[Az, El, Range] = compute_azelrange(r_A_ECEF, satECEF)

%% Problem 3
    % Setup
Code_A = [1, 1, 1, -1, -1, 1, -1];
Code_B = [1, -1, -1, 1, 1, -1, 1];

    % Part a
nShift = 0:length(Code_A);
autoCorrA = [];
autoCorrB = [];
crossCorr = [];
for k = nShift
    autoCorrA = [autoCorrA; calcPRNAutoCorr(Code_A, k)];
    autoCorrB = [autoCorrB; calcPRNAutoCorr(Code_B, k)];
    crossCorr = [crossCorr; calcPRNCrossCorr(Code_A, Code_B, k)];
end

figure
hold on; grid on; grid minor;
plot(nShift, autoCorrA*length(Code_A), 'r-');
plot(nShift, autoCorrB*length(Code_B), 'm-')
plot(nShift, crossCorr*length(Code_A), 'b-');




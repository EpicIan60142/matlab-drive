%% ASEN 5010 Project Task 2 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\..\Utilities\')
addpath('..\Task1\')

R_Mars = 3396.19; % km
h = 400; % km
radius = R_Mars + h; % km

w_0 = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_0 = deg2rad([20; 30; 60]); % Omega, i, theta

x_0 = [radius; EA_0; w_0];

t0 = 0;
dt = 0.5;
tf = 6500;

%% Propagate orbit and find HN
out = RK4_Orbit(x_0, t0, dt, tf);

t = 300;
HN = calcHN(t, out)

%% Save answer to text file
f1 = fopen("HN_ans.txt",'w');
ans_HN = fprintf(f1, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", HN(1,1), HN(1,2), HN(1,3), HN(2,1), HN(2,2), HN(2,3), HN(3,1), HN(3,2), HN(3,3));
fclose(f1);




%% ASEN 5010 Task 4 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\..\Utilities')
addpath('..\Task1')
addpath('..\Task2')

R_Mars = 3396.19; % km
h = 400; % km
radius = R_Mars + h; % km

w_0 = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_0 = deg2rad([20; 30; 60]); % Omega, i, theta

x_0 = [radius; EA_0; w_0];

t0 = 0;
dt = 0.5;
tf = 6500;

%% Propagate orbit and find R_nN
out = RK4_Orbit(x_0, t0, dt, tf);

t = 330;
RnN = calcRnN(t, out)
w_RnN = calcW_RnN(t, out)

%% Save answer to text file
f1 = fopen("RnN_ans.txt",'w');
ans_RnN = fprintf(f1, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", RnN(1,1), RnN(1,2), RnN(1,3), RnN(2,1), RnN(2,2), RnN(2,3), RnN(3,1), RnN(3,2), RnN(3,3));
fclose(f1);

f2 = fopen("w_RnN_ans.txt",'w');
ans_w_RnN = fprintf(f2, "%.8f %.8f %.8f", w_RnN(1), w_RnN(2), w_RnN(3));
fclose(f2);






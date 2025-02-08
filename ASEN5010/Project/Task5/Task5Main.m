%% ASEN 5010 Task 5 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\..\Utilities')
addpath('..\Task1')
addpath('..\Task2')
addpath('..\Task4')

R_Mars = 3396.19; % km
h_LMO = 400; % km
h_GMO = 17028.01; % km
radius_LMO = R_Mars + h_LMO; % km
radius_GMO = R_Mars + h_GMO; % km

w_0_LMO = [0; 0; 0.000884797]; % rad/s, O frame coords
EA_0_LMO = deg2rad([20; 30; 60]); % Omega, i, theta

w_0_GMO = [0; 0; 0.0000709003]; % rad/s, O frame coords
EA_0_GMO = deg2rad([0; 0; 250]); % Omega, i, theta

x_0_LMO = [radius_LMO; EA_0_LMO; w_0_LMO];
x_0_GMO = [radius_GMO; EA_0_GMO; w_0_GMO];

t0 = 0;
dt = 0.5;
tf = 6500;

%% Propagate orbits and find R_cN
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

t = 330;
RcN = calcRcN(t, out_GMO, out_LMO)
w_RcN = calcW_RcN(t, dt, out_GMO, out_LMO)

%% Save answer as text file
f1 = fopen("RcN_ans.txt",'w');
ans_RcN = fprintf(f1, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f", RcN(1,1), RcN(1,2), RcN(1,3), RcN(2,1), RcN(2,2), RcN(2,3), RcN(3,1), RcN(3,2), RcN(3,3));
fclose(f1);

f2 = fopen("w_RcN_ans.txt",'w');
ans_w_RnN = fprintf(f2, "%.8f %.8f %.8f", w_RcN(1), w_RcN(2), w_RcN(3));
fclose(f2);

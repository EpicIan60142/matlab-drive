%% ASEN 5010 Task 6 Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
addpath('..\..\Utilities')
addpath('..\Task1')
addpath('..\Task2')
addpath('..\Task3')
addpath('..\Task4')
addpath('..\Task5')

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

sigBN_0 = [0.3; -0.4; 0.5];
omegBN_0 = deg2rad([1.00; 1.75; -2.20]);

t0 = 0;
dt = 0.5;
tf = 6500;

%% Propagate orbits and find R_cN
out_LMO = RK4_Orbit(x_0_LMO, t0, dt, tf);
out_GMO = RK4_Orbit(x_0_GMO, t0, dt, tf);

t = t0;

RsN = calcRsN();
w_RsN = [0; 0; 0];

RnN = calcRnN(t, out_LMO);
w_RnN = calcW_RnN(t, out_LMO);

RcN = calcRcN(t, out_GMO, out_LMO);
w_RcN = calcW_RcN(t, dt, out_GMO, out_LMO);

[sigBR_sun, omegBR_sun] = calcError(sigBN_0, omegBN_0, RsN, w_RsN) % Sun-pointing error

[sigBR_nadir, omegBR_nadir] = calcError(sigBN_0, omegBN_0, RnN, w_RnN) % Sun-pointing error

[sigBR_GMO, omegBR_GMO] = calcError(sigBN_0, omegBN_0, RcN, w_RcN) % Sun-pointing error

%% Formulate text files for submission
f1 = fopen("sigBR_sun_ans.txt",'w');
ans_sigBR_sun = fprintf(f1, "%.8f %.8f %.8f", sigBR_sun(1), sigBR_sun(2), sigBR_sun(3));
fclose(f1);

f2 = fopen("omegBR_sun_ans.txt",'w');
ans_omegBR_sun = fprintf(f2, "%.8f %.8f %.8f", omegBR_sun(1), omegBR_sun(2), omegBR_sun(3));
fclose(f2);


f3 = fopen("sigBR_nadir_ans.txt",'w');
ans_sigBR_nadir = fprintf(f3, "%.8f %.8f %.8f", sigBR_nadir(1), sigBR_nadir(2), sigBR_nadir(3));
fclose(f3);

f4 = fopen("omegBR_nadir_ans.txt",'w');
ans_omegBR_nadir = fprintf(f4, "%.8f %.8f %.8f", omegBR_nadir(1), omegBR_nadir(2), omegBR_nadir(3));
fclose(f4);


f5 = fopen("sigBR_GMO_ans.txt",'w');
ans_sigBR_GMO = fprintf(f5, "%.8f %.8f %.8f", sigBR_GMO(1), sigBR_GMO(2), sigBR_GMO(3));
fclose(f5);

f6 = fopen("omegBR_GMO_ans.txt",'w');
ans_omegBR_GMO = fprintf(f6, "%.8f %.8f %.8f", omegBR_GMO(1), omegBR_GMO(2), omegBR_GMO(3));
fclose(f6);



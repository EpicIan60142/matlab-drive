%% ASEN 4013 HW 4 Problem 3.33
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Solve
gamma = 1.325;

ratio = 926.99/1475.106;

f = @(M) ((2*(gamma + 1)*M.^2) / (1+gamma*M.^2).^2)*(1 + 0.5*(gamma-1)*M.^2) - ratio;

fsolve(f, 0.5)



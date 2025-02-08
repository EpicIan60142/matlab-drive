%% ECEN 4138 HW 8 Problem 5.20
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Define L(s)

s = tf('s');

L1 = (s+3)/(s*(s+1)*(s^2+4*s+5)); % L(s) for 5.20

L = L1; % Choose root locus

%% Plot root locus
rlocus(L)
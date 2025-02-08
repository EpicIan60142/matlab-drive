%% ECEN 4138 HW 7 Problem 5.7abd
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Define L(s)

s = tf('s');

L1 = (s+3)/(s*(s+10)*(s^2+2*s+2)); % L(s) for 5.7a
L2 = (s+3)/(s^2*(s+10)*(s^2+6*s+25)); % L(s) for 5.7b
L3 = ((s+3)*(s^2+4*s+68))/(s^2*(s+10)*(s^2+4*s+85)); % L(s) for 5.7d

L = L3; % Choose root locus

%% Plot root locus
rlocus(L)

%% ECEN 4138 HW 7 Problem 5.8abd
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Define L(s)

s = tf('s');

L1 = (s+2)/((s+10)*(s^2-1)); % L(s) for 5.8a
L2 = (s+2)/((s*(s+10)*(s^2-1))); % L(s) for 5.8b
L3 = (s^2+2*s+1)/(s*(s+20)^2*(s^2-2*s+2)); % L(s) for 5.8d

L = L3; % Choose root locus

%% Plot root locus
rlocus(L)

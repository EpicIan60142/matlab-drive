%% ECEN 4138 HW 10 problem 6.24
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup

s = tf('s');

G = (0.0001/s)*((s+0.80)/(s+0.01))^2*((s+1)/(s+2))*(1/(s+40))^2

minreal(zpk(G))

%% Plotting
figure
bode(G)

figure
rlocus(G)

figure
nyquist(G)




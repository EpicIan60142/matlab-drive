%% ECEN 4138 HW 10 problem 6.21
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup

s = tf('s');

G = 1/((s+1)*(s^2+2*s+2)); % a

minreal(zpk(G))

%% Plotting
figure
bode(G)

figure
rlocus(G)

figure
nyquist(10.89*G)




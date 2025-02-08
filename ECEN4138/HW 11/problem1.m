%% ECEN 4138 HW 11 Problem 6.25
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

% s = 0.001*1i;

G = 1*(s/10 + 1)/(s*(s-1)*(s/100 + 1));
% G = 1/(s-1)

%% Plotting
figure
bode(G)

figure
nyquist(G)

figure
rlocus(G)

figure
step(feedback(G,1))



%% ECEN 4138 HW 11 Problem 6.26
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

G = 25*(s+1)/(s*(s+2)*(s^2 + 2*s + 16));

%% Plotting
figure
bode(G)

figure
nyquist(G)

figure
rlocus(G)

figure
step(feedback(G,1))

figure
margin(G)



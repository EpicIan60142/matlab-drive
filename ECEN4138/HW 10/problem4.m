%% ECEN 4138 HW 10 problem 6.18
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup

s = tf('s');

G = 1/(s^2 + 4); % a

minreal(zpk(G))

%% Plotting
figure
bode(G)

figure
rlocus(G)

figure
nyquist(G)




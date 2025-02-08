%% ECEN 4138 HW 10 problem 6.19
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup

s = tf('s');

% G = (s+2)/(s+10); % a
G = 1/((s+10)*(s+2)^2); % b
% G = ((s+10)*(s+1))/((s+100)*(s+2)^3); % c

minreal(zpk(G))

%% Plotting
figure
bode(G)

figure
rlocus(G)

figure
nyquist(G)




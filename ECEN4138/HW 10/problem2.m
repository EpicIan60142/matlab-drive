%% ECEN 4138 HW 10 problem 6.22
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup

s = tf('s');

G1 = (s+1)/(s+10); % minimum phase system
G2 = -(s-1)/(s+10); % non-minimum phase system
G3 = (s+1)/(s-10); % unstable system

G = G3;
minreal(zpk(G))

%% Plotting
figure
bode(G)

figure
rlocus(G)

figure
nyquist(G)




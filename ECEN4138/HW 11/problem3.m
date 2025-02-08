%% ECEN 4138 HW 11 Problem 6.45
% Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

G = 50000/(s*(s+1)*(s+50));
C = (1/5)*((0.07179*s + 1)/(0.01616*s + 1))^2;

L = C*G;

%% Plotting
figure
margin(L)

figure
margin(G)

figure
bode(C)

% figure
% nyquist(L)
% 
% figure
% rlocus(L)
% 
% figure
% step(feedback(L,1))



%% ECEN 4138 HW 9 Problem 6.4d
% - Ian Faber

%% Housekeeping
clc; clear; close all;

%% Setup
s = tf('s');

L = ((s+2)*(s+8))/(s*(s+1)*(s+5)*(s+10));

w = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8, 10, 20, 50, 80, 100];

%% Make bode plot
figure
bode(L, w)


